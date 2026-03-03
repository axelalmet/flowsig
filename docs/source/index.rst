Welcome to FlowSig's documentation!
===================================

**FlowSig** A python package to infer communication-driven intercellular flows from single-cell and spatial transcriptomics.

Installation
------------

The easiest way to currently install FlowSig is to create a virtual environment and then install within it.::

.. code-block:: bash

   # Create the virtual environment
   python3 -m venv flowsigenv

   # Activate the virtual environment
   source flowsigenv/bin/activate

**N.B. make sure you're using at least Python 3.11!**

1. **First approach:** pip install.

   .. code-block:: bash
      pip install flowsig


2. **Backup approoach:** install directly from GitHub repository.

   .. code-block:: bash
      # Clone the repository
      git clone https://github.com/axelalmet/flowsig.git
      cd ./flowsig/

      # Install
      pip3 install .

.. toctree::
   :maxdepth: 2
   :caption: Tutorials

   tutorials/index

.. toctree::
   :maxdepth: 3
   :caption: API Reference

   autoapi/flowsig/index

.. toctree::
   :caption: Links

   GitHub <https://github.com/axelalmet/flowsig>