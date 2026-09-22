$ErrorActionPreference = 'Stop'
$taskTools = Join-Path $PSScriptRoot '.tools'
New-Item -ItemType Directory -Path $taskTools -Force | Out-Null
$rArchive = Join-Path $taskTools 'R-4.6.1-win.exe'
$extractArchive = Join-Path $taskTools 'innoextract-1.9-windows.zip'

function Get-VerifiedArchive([string]$Url, [string]$Path, [string]$ExpectedMd5) {
    if (-not (Test-Path -LiteralPath $Path)) {
        Invoke-WebRequest -Uri $Url -OutFile $Path -UseBasicParsing
    }
    $actualMd5 = (Get-FileHash -LiteralPath $Path -Algorithm MD5).Hash.ToLower()
    if ($actualMd5 -ne $ExpectedMd5) { throw "Published checksum mismatch for $Path" }
}

Get-VerifiedArchive 'https://cran.r-project.org/bin/windows/base/R-4.6.1-win.exe' $rArchive '7907f3a20ec8ec88cd0da279024b8e27'
Get-VerifiedArchive 'https://constexpr.org/innoextract/files/innoextract-1.9-windows.zip' $extractArchive '72d0d0dd874b6236eaa44411f4470ee1'
$extractDir = Join-Path $taskTools 'innoextract'
if (-not (Test-Path -LiteralPath (Join-Path $extractDir 'innoextract.exe'))) {
    Expand-Archive -LiteralPath $extractArchive -DestinationPath $extractDir
}
$rscript = Join-Path $taskTools 'R-portable/app/bin/Rscript.exe'
if (-not (Test-Path -LiteralPath $rscript)) {
    & (Join-Path $extractDir 'innoextract.exe') --extract --silent --output-dir (Join-Path $taskTools 'R-portable') $rArchive
    if ($LASTEXITCODE -ne 0) { throw 'Portable R extraction failed' }
}
& $rscript --vanilla (Join-Path $PSScriptRoot 'install_dependencies.R') $PSScriptRoot
if ($LASTEXITCODE -ne 0) { throw 'Local R package installation failed' }
