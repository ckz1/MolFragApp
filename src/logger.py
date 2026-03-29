import sys
import logging

logger = logging.getLogger(__package__)
logger.setLevel(logging.DEBUG)

stdout_handler = logging.StreamHandler(stream=sys.stdout)
stdout_handler.setLevel(logging.INFO)
stdout_handler.setFormatter(logging.Formatter(
    fmt='[%(asctime)s] %(levelname)s: %(message)s',
    # datefmt='%H:%M:%S',
))
logger.addHandler(stdout_handler)