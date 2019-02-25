import itertools as itools
import networkx as nx
from pMuTT.reaction import Reactions

class Network(Reactions):
    """Reaction network. Inherits from :class:`~pMuTT.reaction.Reactions`

    Attributes
    ----------
        reactions : list of :class:`~pMuTT.reaction.Reaction` objects
            Formation reactions for each phase. Reactions should be written
            with consistent reference species to obtain meaningful data.
        graph : networkx.DiGraph object
            Graph representing the reaction network. Nodes correspond to species
            and the edges correspond to reactions.
    """

    def __init__(self, reactions):
        super().__init__(reactions=reactions)
        self.update_network()
    
    def update_network(self, include_TS=True, key='name'):
        """Updates the reaction network
        
        Parameters
        ----------
            include_TS : bool, optional
                Whether transition states should be included. Default is True
            key : str, optional
                Attribute to use as the key in the output dictionary. Default is
                name
        """
        self.network = nx.DiGraph()

        self._add_nodes(include_TS=include_TS, key=key)
        self._add_edges(include_TS=include_TS)
        
    def _add_nodes(self, include_TS=True, key='name'):
        """Adds the species as nodes to the reaction network.
        
        Parameters
        ----------
            include_TS : bool, optional
                Whether transition states should be included. Default is True
            key : str, optional
                Attribute to use as the key in the output dictionary. Default is
                name
        """
        species = self.get_species(include_TS=include_TS, key=key)
        for specie_name, specie in species.items():
            self.network.add_node(specie_name, specie=specie)

    def _add_edges(self, include_TS=True, key='name'):
        """Adds the reactions as edges to the reaction network.
        
        Parameters
        ----------
            include_TS : bool, optional
                Whether transition states should be included. Default is True
            key : str, optional
                Attribute to use as the key in the output dictionary. Default is
                name
        """
        for reaction in self.reactions:
            for reactant, product in itools.product(reaction.reactants, 
                                                    reaction.products):
                self.network.add_edge(reactant[key], 
                                      product[key], 
                                      reaction=reaction,
                                      ts=False)
            if include_TS and reaction.transition_state is not None:
                for reactant, ts in itools.product(reaction.reactants, 
                                                   reaction.transition_state):
                    self.network.add_edge(reactant[key], 
                                          ts[key],
                                          reaction=reaction,
                                          ts=True)
                for ts, product in itools.product(reaction.reactants, 
                                                  reaction.products):
                    self.network.add_edge(ts[key], 
                                          product[key], 
                                          ts=True)
        
    def draw_network(self, layout='kamada_kawai'):
        """Draws the reaction network

        Parameters
        ----------
            layout : str, optional
                Layout to use. See `networkx documentation`_ for supported options.
                Default is 'kamada_kawai'

        .. _`networkx documentation`: https://networkx.github.io/documentation/stable/reference/drawing.html#module-networkx.drawing.layout
        """
        pass