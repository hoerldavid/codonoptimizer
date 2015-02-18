import sys, os

def getCwd():
    '''
    as described in cx_Freeze documentation:
    http://cx-freeze.readthedocs.org/en/latest/faq.html#using-data-files
    '''
    if getattr(sys, 'frozen', False):
        # The application is frozen
        return os.path.dirname(sys.executable)
    else:
        # The application is not frozen
        # Change this bit to match where you store your data files:
        return os.path.dirname(__file__)