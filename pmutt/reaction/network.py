from copy import copy
import itertools as itools
import heapq
import numpy as np
import networkx as nx
import pygal
from scipy import interpolate
from matplotlib import pyplot as plt
from pmutt import _get_specie_kwargs, _force_pass_arguments, _is_iterable
from pmutt.reaction import (Reactions, _write_reaction_state,
                            _parse_reaction_state)


class Network(Reactions):
    """Reaction network. Inherits from :class:`~pmutt.reaction.Reactions`

    Attributes
    ----------
        reactions : list of :class:`~pmutt.reaction.Reaction` objects
            Formation reactions for each phase. Reactions should be written
            with consistent reference species to obtain meaningful data.
        graph : networkx.DiGraph object
            Graph representing the reaction network. Nodes correspond to
            reaction states and the edges correspond to reactions. Nodes are
            identified using a frozen set of tuples where the first element is
            the name of the species and the second element is the stoichiometry.
            e.self.graph. H2 + 0.5O2 is represented as
            frozenset([('H2', 1.), ('O2', 0.5)])

            The nodes have several attributes:
            
            - name
            - species
            - stoich
            - is_transition_state
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
                Attribute to use as the key in the output dictionary.
                Default is name
        """
        self.graph = nx.Graph()
        for reaction in self.reactions:
            states = [reaction.reactants, reaction.products]
            states_stoich = [reaction.reactants_stoich,reaction.products_stoich]

            # Add transition state if available
            if reaction.transition_state is not None:
                states.append(reaction.transition_state)
                states_stoich.append(reaction.transition_state_stoich)

            # Add nodes
            state_sets = []
            for i, (state, stoich) in enumerate(zip(states, states_stoich)):
                species_set = state_to_set(state, stoich)
                state_sets.append(species_set)

                # Transition states occupy index 2 of states
                is_transition_state = (i == 2)
                self.graph.add_node(species_set, species=state, stoich=stoich,
                                    is_transition_state=is_transition_state,
                                    name=_write_reaction_state(species=state,
                                                               stoich=stoich))
            # Add edges
            if reaction.transition_state is None:
                self.graph.add_edge(state_sets[0], state_sets[1])
            else:
                self.graph.add_edge(state_sets[0], state_sets[2])
                self.graph.add_edge(state_sets[1], state_sets[2])

    def get_min_E_span(self, source, target, units=None, species_delimiter='+',
                       cutoff=None, **kwargs):
        # Convert target and source strings to sets
        if source is None:
            source_set = set()
        else:
            source_names, source_stoich = _parse_reaction_state(
                    reaction_str=source, species_delimiter=species_delimiter)
            source_set = state_str_to_set(species_names=source_names,
                                          stoich=source_stoich)

        target_sets = _get_target_sets(target=target,
                                       species_delimiter=species_delimiter)
        E_spans = []
        for path in nx.all_simple_paths(self.graph, source=source_set,
                                        target=target_sets, cutoff=cutoff):
            E_spans.append(self.get_E_span(path=path, units=units, **kwargs))
        return np.amin(E_spans)

    def get_E_span(self, path, units=None, **kwargs):
        """Gets the energy span of a set of reactions. Equations sourced from

        * Kozuch, S.; Shaik, S. How to Conceptualize Catalytic Cycles? The
          Energetic Span Model. Acc. Chem. Res. 2011, 44 (2), 101–110. 
          https://doi.org/10.1021/ar1000956.


        
        :math:`\\delta E = T_{TDTS} - I_{TDI}`
        
        if the TOF-determining transition state (TSTS) appears after the
        TOF-determining intermediate (TDI):

        :math:`\\delta E = T_{TDTS} - I_{TDI} + \\Delta G_r`

        if the TSTS appears before the TDI:

        Parameters
        ----------
            path : list
                Pathway from networkx.get_simple_paths where each element is
                the name of a node along a path.
            units : str
                Units as string. See :func:`~pmutt.constants.R` for accepted
                units but omit the '/K' (e.g. J/mol).
            kwargs : keyword arguments
                Parameters to evaluate Gibbs energy at each state.
        Returns
        -------
            E_span : float
                Energy span of the pathway
        """
        # Get Gibbs energy for each state along path
        G = []
        for state in path:
            species = self.graph.nodes[state]['species']
            stoich = self.graph.nodes[state]['stoich']
            if units is None:
                G.append(get_state_quantity(species=species, stoich=stoich,
                                            method_name='get_GoRT', **kwargs))
            else:
                G.append(get_state_quantity(species=species, stoich=stoich,
                                            method_name='get_G', units=units,
                                            **kwargs))
        # Get indices for TDI and TDTS
        min_i = np.argmin(G)
        max_i = np.argmax(G)

        energy_span = G[max_i] - G[min_i]
        # If the TDTS is before the TDI, add the Gibbs change of the cycle
        if max_i < min_i:
            energy_span += G[-1] - G[0]
        return energy_span

    def plot_network(self, layout='kamada_kawai_layout', source=None,
                     target=None, species_delimiter='+'):
        """Draws the reaction network

        Parameters
        ----------
            layout : str, optional
                Layout to use. See `networkx documentation`_ for supported
                options. Default is 'kamada_kawai_layout'
            source : str, optional
                Initial state as string. This node will be colored red. It not
                specified, no nodes will be colored red
            target : str or list of str, optional
                Final state as string. This node will be colored green. If not
                specified, no nodes will be colored green
            species_delimiter : str, optional
                Delimiter separating species in ``source and ``target``.
                Default is '+'
        Returns
        -------
            figure : `matplotlib.figure.Figure`_
                Add plot to this figure. If not specified, one will be
                generated
            axes : `matplotlib.axes.Axes.axis`_, optional
                Adds plot to this axis. If not specified, one will be generated

        .. _`networkx documentation`: https://networkx.github.io/documentation/stable/reference/drawinself.graph.html#module-networkx.drawinself.graph.layout
        .. _`matplotlib.figure.Figure`: https://matplotlib.org/api/_as_gen/matplotlib.figure.Figure.html
        .. _`matplotlib.axes.Axes.axis`: https://matplotlib.org/api/_as_gen/matplotlib.axes.Axes.axis.html
        """
        colors = []
        labels = {}
        layout = getattr(nx, layout)
        pos = layout(self.graph)


        # Convert target and source strings to sets
        if source is None:
            source_set = set()
        else:
            source_names, source_stoich = _parse_reaction_state(
                    reaction_str=source, species_delimiter=species_delimiter)
            source_set = state_str_to_set(species_names=source_names,
                                          stoich=source_stoich)

        target_sets = _get_target_sets(target=target,
                                       species_delimiter=species_delimiter)

        # Assign colors and labels
        for node in self.graph.nodes(data=True):
            node_set = node[0]
            is_transition_state = node[1]['is_transition_state']
            if node_set == source_set:
                color = 'red'
            elif node_set in target_sets:
                color = 'green'
            elif is_transition_state:
                color = 'skyblue'
            else:
                color = 'gray'
            colors.append(color)
            labels[node[0]] = node[1]['name']
        figure, axes = plt.subplots()
        nx.draw_networkx(self.graph, pos=pos, node_color=colors, ax=axes,
                         with_labels=False)
        nx.draw_networkx_labels(self.graph, pos=pos, labels=labels, ax=axes)
        return figure, axes

    def plot_coordinate_diagram(self, source, target, method_name, units=None,
                                cutoff=None, max_energy_span=None,
                                max_paths=None, pathway_numbers=None,
                                min_x_spacing=1., x_width=0.5,
                                x_scale_TS=0.5, y_scale_TS=0.5,
                                x_label_offset=-0.1, y_label_offset=0.1,
                                species_delimiter='+', viewer='matplotlib',
                                show_state_table=True, show_state_labels=True,
                                table_font_size=None, table_width_ratio=[3, 1],
                                show_energy_span=False,
                                energy_span_format='.2f', colors=None,
                                **kwargs):
        """Plots the reaction coordinate diagram

        Parameters
        ----------
            source : str
                Initial state as string. All pathways will start here.
            target : str or list of str
                Final state as string. All pathways will end here
            method_name : str
                Method to evaluate property of the states. Examples include:
                'get_H', 'get_HoRT', 'get_G', 'get_GoRT'
            units : str, optional
                Units to use to evaluate method_name. Must be specified if the
                `method_name` returns a dimensional property
            cutoff : int, optional
                Maximum number of states in the pathway. If not specified, all
                pathways are shown
            max_energy_span : float, optional
                If specified, pathways with larger energy spans are eliminated
            max_paths : int, optional
                If specified, the number of pathways plotted are limited
            pathway_numbers : list of int, optional
                If specified, only certain pathways are plotted            
            min_x_spacing : float, optional
                Minimum spacing between states. Default is 1.
            x_width : float, optional
                Spacing of stable states. Default is 0.5
            x_scale_TS : float, optional
                Value between 0 and 1 that controls curvature of transition
                state peaks. Higher values produce sharper peaks. Default is
                0.5
            y_scale_TS : float, optional
                Value between 0 and 1 that controls curvature of transition
                state peaks. Higher values produce sharper peaks. Default is
                0.5
            x_label_offset : float, optional
                Horizontal value to offset TS_label from the TS position. This
                value scales with the difference between major ticks. Negative
                values will shift the label leftward. Default is -0.1
            y_label_offset : float, optional
                Vertical value to offset TS_label from the TS position. This
                value scales with the difference between major ticks. Negative
                values will shift the label downwards. Default is 0.1
            species_delimiter : str, optional
                Delimiter that separate species for target and source.
                Leading and trailing spaces will be trimmed. Default is '+'
            viewer : str, optional
                Visualization package to use. Currently, the accepted options
                are: 

                - 'matplotlib' (default)
                - 'pygal'
            show_state_table : bool, optional
                Only applies if `viewer` = 'matplotlib'. If True, a table of
                the states is printed with the diagram. Default is True
            show_state_labels : bool, optional
                Only applies if `viewer` = 'matplotlib'. If True, numbers are
                added to the states which correspond to the entries in the
                table. Default is True
            table_font_size : int, optional
                Only applies if `viewer` = 'matplotlib'. Controls the text font
                size. If not specified, font rescales with figure size
            table_width_ratio : list of int, optional
                Only applies if `viewer` = 'matplotlib'. Controls the relative
                size of diagram to table. i.e. [2, 1] will make the diagram
                width twice as large as the table. Default is [3, 1]
            show_energy_span : bool, optional
                If True, adds energy span value to legend. Default is True
            energy_span_format : str, optional
                String format for energy span in legend. Default is 2 floating
                decimal points
            colors : list of str, optional
                Colors to use for reaction plots
            kwargs: keyword arguments
                Extra arguments that will be fed to evaluate reaction states
        Returns
        -------
            figure : `matplotlib.figure.Figure`_
                Figure
            axes : tuple of `matplotlib.axes.Axes.axis`_
                Axes of the plot.
        Raises
        ------
            ValueError : Raised when `viewer` is not supported.
        """
        # Get the y axis value for the axis label
        y_title = method_name.replace('get_', '')
        if units is not None:
            y_title = '{} ({})'.format(y_title, units)

        # Initialize plot using appropriate viewer
        if viewer == 'matplotlib':
            # Split graph into two axes if including the table
            if show_state_table:
                fig, axes = plt.subplots(
                        ncols=2,
                        gridspec_kw={'width_ratios': table_width_ratio})
            else:
                fig, axes = plt.subplots()
                axes = [axes]
        elif viewer == 'pygal':
            # Use the tooltip x value to indicate what the y value indicates
            x_value_formatter = lambda x: y_title

            # Edit style sheet to have the same color for points and colors
            style = pygal.style.DefaultStyle
            new_colors = []
            if colors is None:
                colors = style.colors
            for color in colors:
                new_colors.extend([color]*2)
            style.colors = new_colors

            # Initialize the graph
            graph = pygal.XY(x_title='Reaction Coordinate',
                             y_title=y_title,
                             pretty_print=True, show_y_guides=False,
                             show_x_guides=False, show_x_labels=False,
                             x_value_formatter=x_value_formatter, style=style,
                             truncate_legend=-1)
        else:
            raise ValueError('Viewer {} not supported. Type '
                             'help(pmutt.reaction.network.Network) '
                             'for supported options.'.format(viewer))

        # If the pathway to plot was specified as an integer, convert to a list
        if pathway_numbers is not None and not _is_iterable(pathway_numbers):
            pathway_numbers = [pathway_numbers]

        # Encode inital and final node
        source_names, source_stoich = _parse_reaction_state(source)
        source_set = state_str_to_set(source_names, source_stoich)
        target_sets = _get_target_sets(target=target,
                                       species_delimiter=species_delimiter)

        # Get all the pathways and associated data for sorting
        paths = list(path for path in nx.all_simple_paths(self.graph,
                                                          source=source_set,
                                                          target=target_sets,
                                                          cutoff=cutoff))
        path_lens = list(len(path) for path in paths)
        energy_spans = list(self.get_E_span(path, units, **kwargs)
                            for path in paths)

        # Get n paths with smallest energy span
        if max_paths is not None:
            paths_data = [(span, path_len, path) for span, path_len, path
                          in zip(energy_spans, path_lens, paths)]
            paths_data_reduced = heapq.nsmallest(max_paths, paths_data)
            energy_spans, path_lens, paths = map(list, zip(*paths_data_reduced))

        # Remove pathways greater than the limit if any
        if max_energy_span is not None:
            for i in range(len(paths)-1, -1, -1):
                if energy_spans[i] > max_energy_span:
                    del paths[i], path_lens[i], energy_spans[i]

        # Sort in descending order of path length
        _, paths_sorted, energy_spans_sorted = zip(*sorted(zip(path_lens,
                                                               paths,
                                                               energy_spans),
                                                   reverse=True))
        
        # Determine x values for each state
        x_vals = {}
        for path in paths_sorted:
            x_spacing = min_x_spacing
            # Find duplicates
            duplicates = tuple(state for state in path if state in x_vals)
            # If there are no duplicates, assign all the states to values and
            # move to next path
            if len(duplicates) == 0:
                for i, state in enumerate(path):
                    x_vals[state] = i*x_spacing
                continue
            # If there is only one duplicate
            if len(duplicates) == 1:
                # If it is the last index, then assign to the length of the
                # reaction.
                if duplicates[0] == path[-1]:
                    x_vals[duplicates[0]] = i*x_spacing*(np.max(path_lens)-1)
                else:
                    i = path.index(duplicates[0])
                    prev_state = path[i-1]
                    next_state = path[i+1]
                    x_vals[duplicates[0]] = np.mean([x_vals[prev_state],
                                                     x_vals[next_state]])
                continue
            # Skip this path if all the states are duplicates
            if len(duplicates) == len(path):
                continue

            # Check adjacent duplicates for intermedient elements
            for duplicate_i, duplicate_j in zip(duplicates, duplicates[1:]):
                i = path.index(duplicate_i)
                j = path.index(duplicate_j)
                # Skip if adjacent duplicates are also adjacent in reaction path
                if j - i == 1:
                    continue
                # Calculate spacing
                x_spacing = (x_vals[duplicate_j] - x_vals[duplicate_i])/(j - i)

                # Assign x positions for new states
                x_initial = x_vals[duplicate_i]
                for l, k in enumerate(range(i+1, j), start=1):
                    x_vals[path[k]] = x_spacing*l + x_initial

        # Sort ascending order by energy span
        energy_spans_sorted, paths_sorted = zip(*sorted(zip(energy_spans_sorted, 
                                                            paths_sorted)))

        # Assign x, y values for plot
        labels_list = []
        labels_set = set()
        y_states = {}
        n_paths = len(paths_sorted)
        for i, (path, energy_span) in enumerate(zip(paths_sorted,
                                                    energy_spans_sorted),
                                                start=1):
            # If pathway_numbers set, skips pathways not specified
            if pathway_numbers is not None and i not in pathway_numbers:
                continue

            # Initialize x, y points for continuous line
            x_plot = []
            y_plot = []
            # Initialize x, y points for interactive points (when viewer is
            # pygal)
            x_points = []
            y_points = []

            # Generate legend for trend
            if show_energy_span:
                if units is None:
                    units_str = ''
                else:
                    units_str = units
                path_name = 'Pathway {:>3} ({:%s} {})'%energy_span_format
                path_name = path_name.format(i, energy_span, units_str)
            else:
                path_name = 'Pathway {:>3}'.format(i)
            point_name = 'Points {:>4}'.format(i)

            for j, state in enumerate(path):
                # If unique state found, add it to the label set
                if state not in labels_set:
                    labels_list.append(state)
                    labels_set.add(state)
                # Get x and y value
                x_state = x_vals[state]
                species = self.graph.nodes[state]['species']
                stoich = self.graph.nodes[state]['stoich']
                y_val = get_state_quantity(species=species, stoich=stoich,
                                           method_name=method_name, units=units,
                                           **kwargs)
                # Subtract the initial state's energy
                if j == 0:
                    y_ref = y_val
                y_state = y_val - y_ref
                y_states[state] = y_state

                # Generate continuous points for plot
                if self.graph.nodes[state]['is_transition_state']:
                    # Calculate product properties for y interpolation
                    products = self.graph.nodes[path[j+1]]['species']
                    prod_stoich = self.graph.nodes[path[j+1]]['stoich']
                    y_prod = get_state_quantity(species=products,
                                                stoich=prod_stoich,
                                                method_name=method_name,
                                                units=units,
                                                **kwargs) - y_ref
                    # Fit spline
                    delta_x = x_state - x_plot[-1]
                    delta_y = y_state - y_plot[-1]
                    x_fit = np.array([x_plot[-1],
                                      x_state - delta_x*x_scale_TS,
                                      x_state,
                                      x_state + delta_x*x_scale_TS,
                                      x_state + delta_x])
                    y_fit = np.array([y_plot[-1],
                                      y_state - delta_y*y_scale_TS,
                                      y_state,
                                      (y_state-y_prod)*y_scale_TS+y_prod,
                                      y_prod])
                    tck = interpolate.splrep(x_fit, y_fit, k=2)
                    # Calculate new x and y points from spline fit
                    x_spline = np.linspace(x_state-delta_x, x_state+delta_x,
                                           100)
                    y_spline = interpolate.splev(x_spline, tck)

                    # Get x value corresponding to peak for pygal
                    max_i = np.argmax(y_spline)
                    x_points.append(x_spline[max_i])
                    y_points.append(y_spline[max_i])
    
                    # Add new data to the appropriate lists
                    x_plot.extend(x_spline)
                    y_plot.extend(y_spline)
                else:
                    # For intermediates, use a straight line
                    x_plot.extend([x_state-x_width/2.,
                                   x_state,
                                   x_state+x_width/2.])
                    y_plot.extend([y_state, y_state, y_state])
                    x_points.append(x_state)
                    y_points.append(y_state)
            # Add data to plot
            if viewer == 'matplotlib':
                axes[0].plot(x_plot, y_plot, label=path_name, zorder=n_paths-i)
            elif viewer == 'pygal':
                # Add line
                line_data = [{'value': (x, y),} for x, y in zip(x_plot, y_plot)]
                graph.add(path_name, line_data, show_dots=False)
                # Add interactive points
                point_data = [{'value': (x, y),
                               'label': self.graph.nodes[state]['name'],}
                              for x, y, state in zip(x_points, y_points, path)]
                graph.add(point_name, point_data, stroke=False)

        if viewer == 'matplotlib':
            # Add other misc labels
            axes[0].legend()
            axes[0].set_ylabel(y_title)
            axes[0].set_xlabel('Reaction coordinate')
            axes[0].tick_params(axis='x', which='both', bottom=False, top=False,
                                labelbottom=False)
            # Add state labels
            if show_state_labels:
                for i, label in enumerate(labels_list, start=1):
                    axes[0].text(x=x_vals[label]+x_label_offset,
                                y=y_states[label]+y_label_offset,
                                s='{:^}'.format(i))
            # Add table
            if show_state_table:
                axes[1].axis('off')
                # Setting up table info
                columns = ('State',)
                rows = range(1, len(labels_list)+1)
                cellText = tuple([self.graph.nodes[state]['name']] 
                                for state in labels_list)
                # Adding table
                table = axes[1].table(cellText=cellText, colLabels=columns,
                                      rowLabels=rows, loc='center')
                # Adjust font size
                if table_font_size is not None:
                    table.auto_set_font_size(False)
                    table.set_fontsize(table_font_size)
            return fig, axes
        else:
            return graph


def state_to_set(species, stoich):
    stoich_out = copy(stoich)
    species_names = []
    for specie in species:
        species_names.append(specie.name)
        try:
            specie.reaction
        except AttributeError:
            pass
        else:
            for state, state_stoich in zip((specie.reaction.reactants,
                                            specie.reaction.products),
                                           (specie.reaction.reactants_stoich,
                                            specie.reaction.products_stoich)):
                species_names.append(state_to_set(state, state_stoich))
                stoich_out.append(1)
    return state_str_to_set(species_names, stoich_out)

def _get_target_sets(target, species_delimiter='+'):
    target_sets = []
    if target is None:
        target_sets.append(set())
    else:
        if not _is_iterable(target):
            target = [target]
        for reaction_str in target:
            target_names, target_stoich = _parse_reaction_state(
                    reaction_str=reaction_str,
                    species_delimiter=species_delimiter)
            target_set = state_str_to_set(
                    species_names=target_names, stoich=target_stoich)
            target_sets.append(target_set)
    return target_sets


def state_str_to_set(species_names, stoich):
    return frozenset([(x, y) for x, y in zip(species_names, stoich)])

def get_state_quantity(species, stoich, method_name, **kwargs):
    if method_name == 'get_q':
        state_quantity = 1.
    else:
        state_quantity = 0.

    for specie, coeff in zip(species, stoich):
        # Process the inputs and methods for each specie
        specie_kwargs = _get_specie_kwargs(specie.name, **kwargs)

        method = getattr(specie, method_name)
        if method_name == 'get_q':
            state_quantity *= \
                    _force_pass_arguments(method, **specie_kwargs)**coeff
        else:
            state_quantity += \
                    _force_pass_arguments(method, **specie_kwargs)*coeff
    return state_quantity
