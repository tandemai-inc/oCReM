import logging


class Logger(object):
    def __init__(self, name):
        self._logger = logger = logging.getLogger(name)
        logger.setLevel(logging.DEBUG)

        # create console handler and set level to debug
        ch = logging.StreamHandler()
        ch.setLevel(logging.DEBUG)

        # create formatter
        formatter = logging.Formatter("%(asctime)s - %(name)s - %(levelname)s - %(message)s")

        # add formatter to ch
        ch.setFormatter(formatter)

        # add ch to logger
        logger.addHandler(ch)
        logger.propagate = False

    def info(self, msg):
        """Log *msg* with severity 'INFO'."""
        self._logger.info(msg)

    def critical(self, msg):
        """Log *msg* with severity 'CRITICAL'."""
        self._logger.critical(msg)

    def debug(self, msg):
        """Log *msg* with severity 'DEBUG'."""
        self._logger.debug(msg)

    def warning(self, msg):
        """Log *msg* with severity 'WARNING'."""
        self._logger.warning(msg)

    def error(self, msg):
        """Log *msg* with severity 'ERROR' and terminate with status 2."""
        self._logger.error(msg)


LOGGER = Logger("tagen")
