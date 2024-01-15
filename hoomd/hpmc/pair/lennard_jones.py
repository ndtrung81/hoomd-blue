# Copyright (c) 2009-2023 The Regents of the University of Michigan.
# Part of HOOMD-blue, released under the BSD 3-Clause License.

"""Lennard Jones pair potential.

.. invisible-code-block: python

    simulation = hoomd.util.make_example_simulation()
    sphere = hoomd.hpmc.integrate.Sphere()
    sphere.shape['A'] = dict(diameter=0.0)
    simulation.operations.integrator = sphere

    pair =  hoomd.hpmc.pair.LennardJones()
    pair.params[('A', 'A')] = dict(epsilon=1, sigma=1, r_cut=2.5)

    logger = hoomd.logging.Logger()
"""

import hoomd

from .pair import Pair


class LennardJones(Pair):
    """Lennard-Jones pair potential (HPMC).

    Args:
        default_r_cut (float): Default cutoff radius :math:`[\\mathrm{length}]`.
        default_r_on (float): Default XPLOR on radius
          :math:`[\\mathrm{length}]`.
        default_mode (str): Default energy
          shifting/smoothing mode.

    `LennardJones` computes the Lennard-Jones pair potential between every pair
    of particles in the simulation state. The functional form of the potential,
    including its behavior under shifting modes, is identical to that in
    the MD pair potential `hoomd.md.pair.LJ`.

    See Also:
        `hoomd.md.pair.LJ`

        `hoomd.md.pair`

    .. rubric:: Example

    .. code-block:: python

        lennard_jones =  hoomd.hpmc.pair.LennardJones()
        lennard_jones.params[('A', 'A')] = dict(epsilon=1, sigma=1, r_cut=2.5)
        simulation.operations.integrator.pair_potentials = [lennard_jones]

    .. py:attribute:: params

        The potential parameters. The dictionary has the following keys:

        * ``epsilon`` (`float`, **required**) -
          Energy well depth :math:`\\varepsilon` :math:`[\\mathrm{energy}]`.
        * ``sigma`` (`float`, **required**) -
          Particle size :math:`\\sigma` :math:`[\\mathrm{length}]`.
        * ``r_cut`` (`float`): Cutoff radius :math:`[\\mathrm{length}]`.
          Defaults to the value given in ``default_r_cut`` on construction.
        * ``r_on`` (`float`): XPLOR on radius :math:`[\\mathrm{length}]`.
          Defaults to the value given in ``default_r_on`` on construction.
        * ``mode`` (`str`): The energy shifting/smoothing mode: Possible values:
          ``"none"``, ``"shift"``, ``"xplor"``. Defaults to the value given
          in ``default_mode`` on construction.

        Type: `TypeParameter` [`tuple` [``particle_type``, ``particle_type``],
        `dict`]
    """
    _cpp_class_name = "PairPotentialLennardJones"

    def __init__(self,
                 default_r_cut=None,
                 default_r_on=0.0,
                 default_mode='none'):
        if default_r_cut is None:
            default_r_cut = float
        else:
            default_r_cut = float(default_r_cut)

        params = hoomd.data.typeparam.TypeParameter(
            'params', 'particle_types',
            hoomd.data.parameterdicts.TypeParameterDict(
                epsilon=float,
                sigma=float,
                r_cut=default_r_cut,
                r_on=float(default_r_on),
                mode=str(default_mode),
                len_keys=2))
        self._add_typeparam(params)
