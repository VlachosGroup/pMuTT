.. _visual:

Visualization
*************

This page contains various functions to visualize functions.

Plot_1D
=======

.. autofunction:: pmutt.plot_1D

Example
-------

Below, we plot the enthalpy, entropy, and Gibbs energy of a water molecule
as a function of temperature at 1 bar and 10 bar.

.. code:: python

    import numpy as np
    from pmutt import plot_1D
    from pmutt.examples import H2O_statmech
    from matplotlib import pyplot as plt
    
    T = np.linspace(300., 500.) # K
    fig1, ax1 = plot_1D(H2O_statmech, x_name='T', x_values=T,
                        methods=('get_H', 'get_S', 'get_G'),
                        get_H_kwargs={'units':'kcal/mol'},
                        get_S_kwargs={'units':'cal/mol/K'},
                        get_G_kwargs={'units': 'kcal/mol'})
    
    # Passing the figure and ax arguments superimpose the plots
    fig1, ax1 = plot_1D(H2O_statmech, x_name='T', x_values=T,
                        P=10., methods=('get_H', 'get_S', 'get_G'),
                        get_H_kwargs={'units':'kcal/mol'},
                        get_S_kwargs={'units':'cal/mol/K'},
                        get_G_kwargs={'units': 'kcal/mol'},
                        figure=fig1, ax=ax1)
    
    # Add legend to matplotlib axes
    ax1[2].legend(['1 bar', '10 bar'])
    ax1[0].set_ylabel('H (kcal/mol)')
    ax1[1].set_ylabel('S (cal/mol/K)')
    ax1[2].set_ylabel('G (kcal/mol)')
    plt.show()

.. image:: ./plot_1D_example.svg

Plot_2D
=======

Below, we plot the enthalpy, entropy, and Gibbs energy of a water molecule
as a function of temperature and pressure.

.. autofunction:: pmutt.plot_2D

Example
-------

.. code:: python

    import numpy as np
    from pmutt import plot_2D
    from pmutt.examples import H2O_statmech
    from matplotlib import pyplot as plt
    
    T = np.linspace(300., 500.) # K
    P = np.logspace(-3, 3) # bar
    
    fig1, ax1, c1, cbar1 = plot_2D(H2O_statmech, 
                                   x1_name='T', x1_values=T,
                                   x2_name='P', x2_values=P,
                                   methods=('get_H', 'get_S', 'get_G'),
                                   get_H_kwargs={'units':'kcal/mol'},
                                   get_S_kwargs={'units':'cal/mol/K'},
                                   get_G_kwargs={'units': 'kcal/mol'})
    plt.show()

.. image:: ./plot_2D_example.png