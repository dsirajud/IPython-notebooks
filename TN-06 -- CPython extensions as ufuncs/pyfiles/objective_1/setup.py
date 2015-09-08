"""
    setup.py file for testmodule.c

    Calling
    $python setup.py build_ext --inplace
    will build the extension library in the current file.
    
    See the distutils section of
    'Extending and Embedding the Python Interpreter'
    at docs.python.org for more information.
"""


from distutils.core import setup, Extension

module1 = Extension('test', sources=['testmodule.c'],
                        include_dirs=['/usr/local/lib'])

setup(name = 'test',
        version='1.0',
        description='This is my test package',
        ext_modules = [module1])
