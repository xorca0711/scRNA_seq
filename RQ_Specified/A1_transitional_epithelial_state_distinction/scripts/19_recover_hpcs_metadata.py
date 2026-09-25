"""Read only remote HDF5 /obs, with hard network and runtime bounds."""
from pathlib import Path
from datetime import datetime, timezone
import hashlib
import io
import json
import time
import urllib.request
import h5py
from anndata.io import read_elem

BASE = Path(__file__).resolve().parents[1]
URL = 'https://ftp.ncbi.nlm.nih.gov/geo/series/GSE277nnn/GSE277777/suppl/GSE277777_combined_data.h5ad'


class BoundedRemote(io.RawIOBase):
    def __init__(self, url):
        self.url = url
        self.position = 0
        self.total = 0
        self.requests = 0
        self.started = time.monotonic()
        self.blocks = {}
        self.block_size = 65536
        self.ranges = []
        req = urllib.request.Request(url, method='HEAD')
        with urllib.request.urlopen(req, timeout=45) as r:
            self.size = int(r.headers['Content-Length'])
            self.etag = r.headers.get('ETag')
            self.modified = r.headers.get('Last-Modified')

    def readable(self):
        return True

    def seekable(self):
        return True

    def tell(self):
        return self.position

    def seek(self, offset, whence=0):
        self.position = offset + (self.position if whence == 1 else self.size if whence == 2 else 0)
        if self.position < 0:
            raise ValueError('Negative seek')
        return self.position

    def readinto(self, target):
        end = min(self.size, self.position+len(target))
        data = bytearray()
        while self.position < end:
            block = self.position // self.block_size
            if block not in self.blocks:
                start = block*self.block_size
                stop = min(self.size, start+self.block_size)-1
                if self.total + stop-start+1 > 64*1024**2 or time.monotonic()-self.started > 900:
                    raise RuntimeError('Frozen metadata-recovery budget exhausted')
                headers = {'Range': f'bytes={start}-{stop}'}
                if self.etag:
                    headers['If-Match'] = self.etag
                req = urllib.request.Request(self.url, headers=headers)
                with urllib.request.urlopen(req, timeout=45) as r:
                    if r.status != 206 or r.headers.get('Content-Range') != f'bytes {start}-{stop}/{self.size}':
                        raise RuntimeError('Server did not honor exact byte range')
                    payload = r.read(stop-start+2)
                if len(payload) != stop-start+1:
                    raise RuntimeError('Unexpected range length')
                self.blocks[block] = payload
                self.total += len(payload)
                self.requests += 1
                self.ranges.append(dict(start=start,end_inclusive=stop,sha256=hashlib.sha256(payload).hexdigest()))
            offset = self.position % self.block_size
            part = self.blocks[block][offset:offset+end-self.position]
            data.extend(part)
            self.position += len(part)
        target[:len(data)] = data
        return len(data)


def main():
    output = BASE/'cache/followup_sources/hpcs_combined_obs.csv.gz'
    report = BASE/'reports/hpcs_metadata_recovery.json'
    if output.exists() or report.exists():
        raise SystemExit('Refusing to overwrite metadata recovery')
    record = dict(utc=datetime.now(timezone.utc).isoformat(),url=URL,
        code_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        contract_sha256=hashlib.sha256((BASE/'config/hpcs_metadata_recovery.md').read_bytes()).hexdigest())
    remote = None
    try:
        remote = BoundedRemote(URL)
        with h5py.File(remote, 'r') as h:
            obs = read_elem(h['obs'])
        obs.to_csv(output, compression='gzip')
        record.update(status='completed',rows=len(obs),columns=list(obs.columns),
            output_path=output.relative_to(BASE).as_posix(),output_sha256=hashlib.sha256(output.read_bytes()).hexdigest())
        print('Metadata recovered:',len(obs),'rows; columns:',list(obs.columns))
    except Exception as exc:
        record.update(status='held',error=repr(exc))
        print('Recovery held:',repr(exc))
    finally:
        if remote:
            record.update(remote_bytes=remote.size,etag=remote.etag,last_modified=remote.modified,
                transferred_bytes=remote.total,range_requests=remote.requests,ranges=remote.ranges,
                elapsed_seconds=round(time.monotonic()-remote.started,2))
        report.write_text(json.dumps(record,indent=2)+'\n')


if __name__ == '__main__':
    main()
