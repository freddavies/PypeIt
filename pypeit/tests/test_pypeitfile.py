""" Tests on PypeItFile """

from pathlib import Path

from astropy.table import Table
from IPython import embed

from pypeit import dataPaths
from pypeit.inputfiles import PypeItFile
from pypeit.tests.tstutils import data_output_path


def pieces_of_pypeitfile():
    """
    Bits needed to generate a PypeIt file
    """
    confdict = {'rdx': {'spectrograph': 'keck_hires'}}

    data = Table()
    data['filename'] = ['b1.fits.gz', 'b27.fits.gz']
    data['frametype'] = ['arc', 'science']
    data['exptime'] = [1., 10.]

    file_paths = [data_output_path('')]

    setup_dict = {'Setup A': ' '}

    # Return
    return confdict, data, file_paths, setup_dict


def test_instantiate():
    # Test of instantiation
    confdict, data, file_paths, setup_dict = pieces_of_pypeitfile()
    pypeItFile = PypeItFile(confdict, file_paths, data, setup_dict)
    # Data files                                    
    filenames = pypeItFile.filenames
    assert 'b1.fits.gz' in filenames[0]

    # Frame types
    frame_type_dict = pypeItFile.frametypes
    assert frame_type_dict['b1.fits.gz'] == 'arc'

    # More tests
    assert pypeItFile.setup_name == 'A'


def test_read_pypeit_file():
    # Read the PypeIt file (backwards compatibility)
    ifile = dataPaths.tests.get_file_path('example_pypeit_file.pypeit')
    pypeItFile = PypeItFile.from_file(ifile)
    assert isinstance(pypeItFile.config, dict)


def test_read_backwards_pypeit_file():
    # Read the PypeIt file (backwards compatibility)
    ifile = dataPaths.tests.get_file_path('example_pypeit_file_backwards.pypeit')
    pypeItFile = PypeItFile.from_file(ifile)
    assert isinstance(pypeItFile.config, dict)


def test_write_pypeit_file():
    # Test writing a PypeIt file
    outfile = Path(data_output_path('tmp_file.pypeit')).absolute()
    if outfile.is_file():
        outfile.unlink()

    # Instantiate
    confdict, data, file_paths, setup_dict = pieces_of_pypeitfile()
    pypeItFile = PypeItFile(confdict, file_paths, data, setup_dict)
    # Write
    pypeItFile.write(outfile)
    assert outfile.is_file(), 'PypeIt file not written'

    # Let's read it too
    pypeItFile2 = PypeItFile.from_file(outfile)
    assert dict(pypeItFile2.config) == confdict, 'Config dict not the same after read/write'
    assert all(pypeItFile2.data == data), 'Data table not the same after read/write'
    assert str(pypeItFile2.file_paths[0]) == file_paths[0], \
        'File paths not the same after read/write'
    assert pypeItFile2.setup == setup_dict, 'Setup dict not the same after read/write'

    # Clean up
    outfile.unlink()

