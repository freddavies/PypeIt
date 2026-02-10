"""
Tests on io module
"""
from pathlib import Path

from astropy.table import Table
from IPython import embed

from pypeit import dataPaths
from pypeit import inputfiles
from pypeit.tests import tstutils


def _pypeitfile_components():
    """
    Bits needed to generate a PypeIt file
    """
    confdict = {'rdx': {'spectrograph': 'keck_hires'}}

    data = Table()
    data['filename'] = ['b1.fits.gz', 'b27.fits.gz']
    data['frametype'] = ['arc', 'science']
    data['exptime'] = [1., 10.]

    file_paths = [tstutils.data_output_path('')]

    setup_dict = {'Setup A': ' '}

    # Return
    return confdict, data, file_paths, setup_dict


def test_grab_rawfiles():

    tst_file = Path(tstutils.data_output_path('test.rawfiles')).absolute()
    if tst_file.exists():
        tst_file.unlink()

    # Download and move all the b*fits.gz files into the local package
    # installation
    tstutils.install_shane_kast_blue_raw_data()

    root = Path(tstutils.data_output_path('')).absolute()
    raw_files = [root / 'b11.fits.gz', root / 'b12.fits.gz']
    assert all([f.exists() for f in raw_files]), 'Files missing'

    tbl = Table()
    tbl['filename'] = [r.name for r in raw_files]
    inputfiles.RawFiles(file_paths=[str(root)], data_table=tbl).write(tst_file)

    _raw_files = inputfiles.grab_rawfiles(file_of_files=str(tst_file))
    assert all(f == Path(_f) for f, _f in zip(raw_files, _raw_files)), 'Bad file_of_files read'

    _raw_files = inputfiles.grab_rawfiles(list_of_files=tbl['filename'], raw_paths=[str(root)])
    assert all(f == Path(_f) for f, _f in zip(raw_files, _raw_files)), 'Bad list_of_files read'

    _raw_files = inputfiles.grab_rawfiles(raw_paths=[str(root)], extension='.fits.gz')
    assert len(_raw_files) == 9, 'Found the wrong number of files'
    assert all([str(root / f) in _raw_files for f in tbl['filename']]), 'Missing expected files'

    tst_file.unlink()


def test_instantiate_pypeitfile():
    # Test of instantiation
    confdict, data, file_paths, setup_dict = _pypeitfile_components()
    pypeItFile = inputfiles.PypeItFile(confdict, file_paths, data, setup_dict)
    # Data files                                    
    filenames = pypeItFile.filenames
    assert 'b1.fits.gz' in filenames[0]

    # Frame types
    frame_type_dict = pypeItFile.frametypes
    assert frame_type_dict['b1.fits.gz'] == 'arc'

    # More tests
    assert pypeItFile.setup_name == 'A'


def test_read_pypeitfile():
    # Read the PypeIt file (backwards compatibility)
    ifile = dataPaths.tests.get_file_path('example_pypeit_file.pypeit')
    pypeItFile = inputfiles.PypeItFile.from_file(ifile)
    assert isinstance(pypeItFile.config, dict)


def test_read_backwards_pypeitfile():
    # Read the PypeIt file (backwards compatibility)
    ifile = dataPaths.tests.get_file_path('example_pypeit_file_backwards.pypeit')
    pypeItFile = inputfiles.PypeItFile.from_file(ifile)
    assert isinstance(pypeItFile.config, dict)


def test_write_pypeitfile():
    # Test writing a PypeIt file
    outfile = Path(tstutils.data_output_path('tmp_file.pypeit')).absolute()
    if outfile.is_file():
        outfile.unlink()

    # Instantiate
    confdict, data, file_paths, setup_dict = _pypeitfile_components()
    pypeItFile = inputfiles.PypeItFile(confdict, file_paths, data, setup_dict)
    # Write
    pypeItFile.write(outfile)
    assert outfile.is_file(), 'PypeIt file not written'

    # Let's read it too
    pypeItFile2 = inputfiles.PypeItFile.from_file(outfile)
    assert dict(pypeItFile2.config) == confdict, 'Config dict not the same after read/write'
    assert all(pypeItFile2.data == data), 'Data table not the same after read/write'
    assert str(pypeItFile2.file_paths[0]) == file_paths[0], \
        'File paths not the same after read/write'
    assert pypeItFile2.setup == setup_dict, 'Setup dict not the same after read/write'

    # Clean up
    outfile.unlink()




