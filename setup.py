from setuptools import setup

# Next 3 lines might be needed due to bug https://github.com/pypa/pip/issues/7953 in setuptools <62.0
import site
import sys
site.ENABLE_USER_SITE = "--user" in sys.argv[1:]

if __name__ == "__main__":
    setup()
