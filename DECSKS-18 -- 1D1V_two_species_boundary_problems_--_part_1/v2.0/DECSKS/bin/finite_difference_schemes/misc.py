def HighPrecisionE(number):
    """Converts a number into a string object
    while retaining a chosen degree of precision. This
    is designed to evade the truncation that is involved
    with str() so that outputs can store numbers with high
    precision

    inputs:
    number -- (number)

    outputs:
    string object with chosen precision in scientific notation
    """

    return "%.22e" % number


def strings_for_writeout():
    """creates a dictionary of the strings below used in
    the outfile writing decoration

    inputs:
    None

    outputs:
    table_strings -- (dict) decorating strings used in output writing
    """
    newline = '\n'
    divider = '-------------------------'
    spc = ' '
    colon = ':'
    equaldivider = '========================='
    outfilename_suffix = '_FD_coefficients.dat'
    number_of_lines_to_read = 'number of lines to read: '

    table_strings = dict(newline = newline,
                         divider = divider,
                         spc = spc,
                         colon = colon,
                         equaldivider = equaldivider,
                         outfilename_suffix = outfilename_suffix,
                         number_of_lines_to_read = number_of_lines_to_read
                        )

    return table_strings

def write_header(dn, LTE):

    strings = strings_for_writeout()

    # locally store names for brevity
    newline = strings['newline']
    spc = strings['spc']
    colon = strings['colon']
    equaldivider = strings['equaldivider']
    outfilename_suffix = strings['outfilename_suffix']
    underscore = '_'
    # open file for decorated table
    decorated_out = open('Table_of_finite_difference_schemes_at_const_LTE_decorated.dat','w')

    # header for decorated table output
    derivative = 'f^(' + str(dn) + ')'
    decorated_out.write(derivative + colon + spc + 'LTE = O(z^' + str(LTE)
                          + ')' + '\n')
    decorated_out.write(equaldivider + newline + newline)

    # filename construction for outfile for each derivative dn at chosen LTE
    dn_outname_prefix = 'f' + str(dn)
    dn_outname_LTE = 'LTE' + str(LTE)
    dn_outname = dn_outname_prefix + underscore + dn_outname_LTE + outfilename_suffix
    dn_schemes_out = open(dn_outname, 'w')

    return decorated_out, dn_schemes_out

def write_scheme(outfile, _w, _stencil, high_precision = 'no'):
    """For a provided outfile by the user, write
    the scheme (no plural!) according to the provided finite difference
    weights (_w) and stencil (_stencil)

    inputs:
    outfile -- (file) an opened ('w' mode) file to be written to
    _w -- (list) weights of the finite difference scheme
    _stencil -- (list) corresponding stencil for the scheme

    outputs:
    outfile -- (file) return the same file after writing is complete
    """
    stencil_size = len(_stencil)

    high_precision = high_precision.lower()

    if high_precision == 'yes':
        for i in range(stencil_size):
            outfile.write( HighPrecisionE(_w[i]) + ', '
                           + str(_stencil[i]) + '\n')

    elif high_precision == 'no':
        for i in range(stencil_size):
            outfile.write( '%2.5f' % _w[i] + ', '
                           + '%2d' % _stencil[i] + '\n')

    return outfile

def write_rest(
        decorated_out,
        dn_schemes_out,
        label,
        p,
        dn,
        stencil,
        w):

    strings = strings_for_writeout()
    colon = strings['colon']
    divider = strings['divider']
    newline = strings['newline']
    spc = strings['spc']
    derivative = 'f^(' + str(dn) + ')'
    number_of_lines_to_read = strings['number_of_lines_to_read']

    decorated_out.write(derivative + colon + spc + label + '\n')
    decorated_out.write(divider + newline)
    decorated_out.write('(weights, stencil)' + newline)
    decorated_out = write_scheme(decorated_out, w, stencil)
    decorated_out.write(divider + newline)

    # write scheme to reader function friendly table
    dn_schemes_out.write(label + newline)
    dn_schemes_out.write(number_of_lines_to_read + str(p+dn) + newline)
    dn_schemes_out = write_scheme(dn_schemes_out, w, stencil, high_precision = 'yes')

    return decorated_out, dn_schemes_out

def write_footer(decorated_out):

    strings = strings_for_writeout()
    equaldivider = strings['equaldivider']
    newline = strings['newline']
    decorated_out.write(equaldivider + newline)

    return decorated_out
