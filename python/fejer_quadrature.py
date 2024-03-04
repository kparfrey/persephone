# From https://chaospy.readthedocs.io/en/master/_modules/chaospy/quadrature/fejer_1.html
import numpy

def fejer_1(order):
    """Backend for Fejer type I quadrature."""
    order = int(order)
    if order == 0:
        return numpy.array([0.5]), numpy.array([1.0])
    order += 1

    abscissas = -0.5 * numpy.cos(numpy.pi * (numpy.arange(order) + 0.5) / order) + 0.5

    steps = numpy.arange(1, order, 2)
    length = len(steps)
    remains = order - length

    kappa = numpy.arange(remains)
    beta = numpy.hstack(
        [
            2 * numpy.exp(1j * numpy.pi * kappa / order) / (1 - 4 * kappa**2),
            numpy.zeros(length + 1),
        ]
    )
    beta = beta[:-1] + numpy.conjugate(beta[:0:-1])

    weights = numpy.fft.ifft(beta)
    assert max(weights.imag) < 1e-15
    weights = weights.real / 2.0

    #return abscissas, weights
    # Only return the weights, since the point locations are already set
    return weights 
