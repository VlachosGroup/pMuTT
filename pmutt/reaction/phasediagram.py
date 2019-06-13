# -*- coding: utf-8 -*-
import numpy as np
import matplotlib
from matplotlib import pyplot as plt
from pmutt import constants as c
from pmutt.reaction import Reactions
from pmutt.io.json import json_to_pmutt, remove_class


class PhaseDiagram(Reactions):
    """Generate phase diagrams based on reactions specified. Inherits from
    :class:`~pmutt.reaction.Reactions`

    Attributes
    ----------
        reactions : list of :class:`~pmutt.reaction.Reaction` objects
            Formation reactions for each phase. Reactions should be written
            with consistent reference species to obtain meaningful data.
        norm_factors : (N,) `numpy.ndarray`_ of float, optional
            Used for normalizing Gibbs energies. These factors could be
            surface areas when calculating surface energies or if the
            reactions stoichiometry is not consistent. Default is an array of
            1. It should have the same length as reactions.

    .. _`numpy.ndarray`: https://docs.scipy.org/doc/numpy/reference/generated/numpy.ndarray.html
    """

    def __init__(self, reactions, norm_factors=None):
        super().__init__(reactions=reactions)
        if norm_factors is None:
            self.norm_factors = np.ones(len(reactions))
        else:
            self.norm_factors = norm_factors

    def to_dict(self):
        """Represents object as dictionary with JSON-accepted datatypes

        Returns
        -------
            obj_dict : dict
        """
        obj_dict = super().to_dict()
        obj_dict['class'] = str(self.__class__),
        return obj_dict

    @classmethod
    def from_dict(cls, json_obj):
        """Recreate an object from the JSON representation.

        Parameters
        ----------
            json_obj : dict
                JSON representation
        Returns
        -------
            PhaseDiagram : PhaseDiagram object
        """
        json_obj = remove_class(json_obj)
        json_obj['reactions'] = [json_to_pmutt(reaction)
                                 for reaction in json_obj['reactions']]
        return cls(**json_obj)

    def get_GoRT_1D(self, x_name, x_values, G_units=None, **kwargs):
        """Calculates the Gibbs free energy for all the reactions for 1 varying
        parameter

        Parameters
        ----------
            x_name : str
                Name of variable to vary
            x_values : iterable object
                x values to use
            G_units : str, optional
                Units for G. If None, uses GoRT. Default is None
            kwargs : keyword arguments
                Other variables to use in the calculation
        Returns
        -------
            GoRT : (M, N) `numpy.ndarray`_ of float
                GoRT values. The first index corresponds to the number of
                reactions. The second index corresponds to the conditions
                specified by x_values.
            stable_phases : (N,) `numpy.ndarray`_ of int
                Each element of the array corresponds to the index of the most
                stable phase at the x_values.
        """
        GoRT = np.zeros(shape=(len(self.reactions), len(x_values)))
        for i, (reaction, norm_factor) in enumerate(zip(self.reactions,
                                                        self.norm_factors)):
            for j, x in enumerate(x_values):
                kwargs[x_name] = x
                GoRT[i, j] = reaction.get_delta_GoRT(**kwargs)/norm_factor

                # Add unit corrections
                if G_units is not None:
                    GoRT[i, j] *= c.R('{}/K'.format(G_units))*kwargs['T']
        stable_phases = np.nanargmin(GoRT, axis=1)
        return (GoRT, stable_phases)

    def plot_1D(self, x_name, x_values, G_units=None, **kwargs):
        """Make a 1D phase diagram.

        Parameters
        ----------
            x_name : str
                Name of variable to vary
            x_values : iterable object
                x values to use
            G_units : str, optional
                Units for G. If None, uses GoRT. Default is None
            kwargs : keyword arguments
                Other variables to use in the calculation
        Returns
        -------
            figure : `matplotlib.figure.Figure`_
                Figure
            ax : `matplotlib.axes.Axes.axis`_
                Axes of the plots.

        .. _`matplotlib.figure.Figure`: https://matplotlib.org/api/_as_gen/matplotlib.figure.Figure.html
        .. _`matplotlib.axes.Axes.axis`: https://matplotlib.org/api/_as_gen/matplotlib.axes.Axes.axis.html
        """
        fig, ax = plt.subplots()
        GoRT, stable_phases = self.get_GoRT_1D(x_name=x_name,
                                               x_values=x_values,
                                               G_units=G_units, **kwargs)
        for GoRT_rxn, rxn in zip(GoRT, self.reactions):
            plt.plot(x_values, GoRT_rxn, label=rxn.to_string())
        ax.legend()
        ax.set_xlabel(x_name)
        if G_units is None:
            ax.set_ylabel('G/RT')
        else:
            ax.set_ylabel('G ({})'.format(G_units))
        return (fig, ax)

    def get_GoRT_2D(self, x1_name, x1_values, x2_name, x2_values,
                    G_units=None, **kwargs):
        """Calculates the Gibbs free energy for all the reactions for two
        varying parameters

        Parameters
        ----------
            x1_name : str
                Name of first variable to vary
            x1_values : iterable object
                x1 values to use
            x2_name : str
                Name of second variable to vary
            x2_values : iterable object
                x2 values to use
            G_units : str, optional
                Units for G. If None, uses GoRT. Default is None
            kwargs : keyword arguments
                Other variables to use in the calculation
        Returns
        -------
            GoRT : (M, N, O) `numpy.ndarray`_ of float
                GoRT values. The first index corresponds to the number of
                reactions. The second index corresponds to the conditions
                specified by x_values.
            stable_phases : (N, O) `numpy.ndarray`_ of int
                Each element of the array corresponds to the index of the most
                stable phase at the x_values.
        """
        GoRT = np.zeros(
                shape=(len(self.reactions), len(x1_values), len(x2_values)))
        for i, (reaction, norm_factor) in enumerate(zip(self.reactions,
                                                        self.norm_factors)):
            for j, x1 in enumerate(x1_values):
                kwargs[x1_name] = x1
                for k, x2 in enumerate(x2_values):
                    kwargs[x2_name] = x2
                    GoRT[i, j, k] = \
                        reaction.get_delta_GoRT(**kwargs)/norm_factor
                    # Add unit corrections
                    if G_units is not None:
                        GoRT[i, j, k] *= c.R('{}/K'.format(G_units)) *\
                                             kwargs['T']
        # Take a transpose
        GoRT_T = GoRT.transpose((1, 2, 0))
        stable_phases = np.zeros((len(x1_values), len(x2_values)))
        for i, GoRT_row in enumerate(GoRT_T):
            stable_phases[i, :] = np.nanargmin(GoRT_row, axis=1)

        return GoRT, stable_phases

    def plot_2D(self, x1_name, x1_values, x2_name, x2_values, G_units=None,
                **kwargs):
        """Make a 2D phase diagram.

        Parameters
        ----------
            x1_name : str
                Name of first variable to vary
            x1_values : iterable object
                x1 values to use
            x2_name : str
                Name of first variable to vary
            x2_values : iterable object
                x2 values to use
            G_units : str, optional
                Units for G. If None, uses GoRT. Default is None
            kwargs : keyword arguments
                Other variables to use in the calculation
        Returns
        -------
            figure : `matplotlib.figure.Figure`_
                Figure
            ax : `matplotlib.axes.Axes.axis`_
                Axes of the plots.
            c : `matplotlib.collections.QuadMesh`_
                Heatmap plot
            cbar : `matplotlib.colorbar.Colorbar`_
                Colorbar for plot

        .. _`matplotlib.figure.Figure`: https://matplotlib.org/api/_as_gen/matplotlib.figure.Figure.html
        .. _`matplotlib.axes.Axes.axis`: https://matplotlib.org/api/_as_gen/matplotlib.axes.Axes.axis.html
        .. _`matplotlib.collections.QuadMesh`: https://matplotlib.org/api/collections_api.html#matplotlib.collections.QuadMesh
        .. _`matplotlib.colorbar.Colorbar`: https://matplotlib.org/api/_as_gen/matplotlib.pyplot.colorbar.html
        """
        # Process input data
        x2_mesh, x1_mesh = np.meshgrid(x2_values, x1_values)
        GoRT, stable_phases = self.get_GoRT_2D(x1_name=x1_name,
                                               x1_values=x1_values,
                                               x2_name=x2_name,
                                               x2_values=x2_values,
                                               G_units=G_units, **kwargs)

        fig, ax = plt.subplots()
        # Choosing color palette
        cmap = plt.get_cmap('viridis')
        norm = matplotlib.colors.BoundaryNorm(np.arange(len(self.reactions)+1),
                                              cmap.N)
        # Create colormap
        c = plt.pcolormesh(x1_mesh, x2_mesh, stable_phases, cmap=cmap,
                           norm=norm, vmin=0, vmax=len(self.reactions))
        # Set colorbar
        cbar = fig.colorbar(c, ticks=np.arange(len(self.reactions))+0.5)
        cbar.ax.set_yticklabels(
                [reaction.to_string() for reaction in self.reactions])
        # Set axis labels
        ax.set_xlabel(x1_name)
        ax.set_ylabel(x2_name)
        return (fig, ax, c, cbar)
