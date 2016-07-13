def f(a, **kwargs):

    b = a + 2

    print "b = %g" % b

    if 'vD' in kwargs.keys():
        vD = kwargs['vD']
        assert ( type(vD) != str)
        print "it's true!"
        print vD

    return None
