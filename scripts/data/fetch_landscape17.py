"""Download and verify Landscape17, then extract the malonaldehyde benchmark."""
import hashlib,json,shutil,urllib.request,zipfile
from pathlib import Path

ROOT=Path(__file__).resolve().parents[2]
FOLDER=ROOT/'data/raw/landscape17'
ARCHIVE=FOLDER/'Landscape17.zip'
URL='https://ndownloader.figshare.com/files/57341180'
SIZE=53019632
MD5='322e3d9d9851c469d0f6069ec5fd846d'


def verify(path):
    return path.stat().st_size==SIZE and hashlib.md5(path.read_bytes()).hexdigest()==MD5


def main():
    FOLDER.mkdir(parents=True,exist_ok=True)
    if not ARCHIVE.exists():
        partial=ARCHIVE.with_suffix('.zip.part')
        with urllib.request.urlopen(URL) as source,open(partial,'wb') as destination:
            shutil.copyfileobj(source,destination)
        if not verify(partial):raise ValueError('Landscape17 download failed size/MD5 verification')
        partial.replace(ARCHIVE)
    if not verify(ARCHIVE):raise ValueError('Existing Landscape17 archive failed size/MD5 verification')
    prefix='Landscape17/malonaldehyde/';destination=FOLDER/'extracted';count=0
    with zipfile.ZipFile(ARCHIVE) as archive:
        for member in archive.infolist():
            if not member.filename.startswith(prefix):continue
            target=(destination/member.filename).resolve()
            if destination.resolve() not in target.parents and target!=destination.resolve():
                raise ValueError('Unsafe archive member path')
            archive.extract(member,destination);count+=1
        if archive.testzip() is not None:raise ValueError('Landscape17 ZIP CRC check failed')
    receipt=dict(dataset='Landscape17',version=1,doi='10.6084/m9.figshare.29949230.v1',
        license='CC BY 4.0',download_url=URL,file_id=57341180,size=SIZE,md5=MD5,
        sha256=hashlib.sha256(ARCHIVE.read_bytes()).hexdigest(),archive_crc_passed=True,
        extracted_subset='malonaldehyde',extracted_members=count)
    (FOLDER/'receipt.json').write_text(json.dumps(receipt,indent=2),encoding='utf-8')
    (ROOT/'reports/aimnet2025_transitionnet_v9/data_receipt.json').write_text(
        json.dumps(receipt,indent=2),encoding='utf-8')
    print(json.dumps(receipt,indent=2))


if __name__=='__main__':main()
