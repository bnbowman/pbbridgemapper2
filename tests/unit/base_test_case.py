import os
import sys
import logging
import unittest
import ConfigParser
import shutil
import traceback

from pbcore.util.Process import backticks


log = logging.getLogger(__name__)

_DIR_NAME = os.path.dirname(os.path.abspath(__file__))
_NOSE_REPORT_CFG = os.path.join(_DIR_NAME, '..', 'nose.cfg')

if os.path.exists(_NOSE_REPORT_CFG):
    _NOSE_REPORT_CFG = os.path.abspath(_NOSE_REPORT_CFG)
    print "loading data config from {0}".format(_NOSE_REPORT_CFG)
else:
    raise IOError("Unable to find test CFG file {f}.".format(f=_NOSE_REPORT_CFG))


# Need to rethink this? Not sure we'll have the correct privileges?
outdir = os.path.join(_DIR_NAME, 'output')


def _get_root_data_dir():
    """Get the root level directory which contains all the unittests files"""
    log.info("Loading data config from {f}".format(f=_NOSE_REPORT_CFG))
    section = 'data'
    option = 'root'
    nosecfg = ConfigParser.SafeConfigParser()
    nosecfg.readfp(open(_NOSE_REPORT_CFG), 'r')
    if nosecfg.has_section(section):
        d = nosecfg.get(section, option)
        data_dir = os.path.abspath(d)
        msg = "Loading data from {d}".format(d=data_dir)
        print msg
        log.info(msg)
        return data_dir
    else:
        msg = "Unable to find section:{s} or option:{o} in {f}".format(s=section, o=option, f=nosecfg)
        log.error(msg)
        KeyError(msg)

ROOT_DATA_DIR = _get_root_data_dir()


def run_backticks(cmd):
    """Util to run an integration test via subprocess.

    Run a cmd in a subprocess and return an exit code (int).
    """
    log.info("Running cmd {c}".format(c=cmd))

    output, rcode, emsg = backticks(cmd)

    if rcode != 0:
        log.error(output)
        log.error(emsg)

    return rcode


class LogManager(object):
    _loggingSetup = False

    @staticmethod
    def _isLoggingSetup():
        return LogManager._loggingSetup

    @staticmethod
    def _setLoggingSetup():
        LogManager._loggingSetup = True


class BaseTestCase(unittest.TestCase):

    def setUp(self):
        """
        Before *every* test create, and init the db
        """
        unittest.TestCase.setUp(self)
        if not os.path.exists(outdir):
            os.mkdir(outdir)

    def tearDown(self):
        if os.path.exists(outdir):
            shutil.rmtree(outdir)

    def get_output_dir(self):
        return outdir

    # ------- @Class methods -----------------#
    @classmethod
    def setUpClass(cls):

        try:
            cls.nosecfg = ConfigParser.SafeConfigParser()
            cls.nosecfg.readfp(open(_NOSE_REPORT_CFG), 'r')
            cls.setUpLog()
            cls._check_data_dir()
        except:
            tb = traceback.format_exc()
            log.error(tb)
            print(tb)
            raise

    @classmethod
    def _check_data_dir(cls):
        d = cls.get_data_dir()
        if d is None:
            raise KeyError("No test data root has been defined in nose.cfg.")
        if not os.path.exists(d):
            raise IOError("Test data root does not exist: {d}".format(d=d))

    @classmethod
    def get_data_dir(cls):
        d = cls.get_cfg_opt('root', section='data')
        if d is not None:
            return os.path.expanduser(d)
        return None

    @classmethod
    def setUpLog(cls):
        """
        If a log file is specified in nose.cfg, log to that file. Else,
        log to the console
        """

        if LogManager._isLoggingSetup():
            # Only want to do this once. setUpLog is called betwixt tests
            log.debug('Logging is already setup.')
            return

        if 'debug' == cls.get_cfg_opt('level', section='log', default='debug'):
            level = logging.DEBUG
        else:
            level = logging.INFO

        str_formatter = '[%(levelname)s] %(asctime)-15s [%(name)s %(funcName)s %(lineno)d] %(message)s'
        formatter = logging.Formatter(str_formatter)
        log_file = cls.get_cfg_opt('file', section='log')
        print "log file", log_file
        # import ipdb; ipdb.set_trace()
        if log_file is not None:
            log_file = os.path.expanduser(log_file)
            log_file = os.path.abspath(log_file)
            if os.path.exists(log_file):
                # delete the last log
                os.remove(log_file)
            else:
                # make the parent dirs if they don't exist
                d = os.path.dirname(log_file)
                if not os.path.exists(d):
                    os.makedirs(d)

            handler = logging.FileHandler(log_file)
            handler.setFormatter(formatter)
            handler.setLevel(level)
            log.addHandler(handler)
        else:
            handler = logging.StreamHandler(sys.stdout)
            handler.setFormatter(formatter)
            handler.setLevel(level)
            log.addHandler(handler)

        log.debug('Set up logger')

        LogManager._setLoggingSetup()

    @classmethod
    def get_cfg_opt(cls, option, section='DEFAULT', default=None):
        '''Get a config option. Returns default, if no such option exists'''
        try:
            return cls.nosecfg.get(section, option)
        except ConfigParser.NoOptionError:
            return default
