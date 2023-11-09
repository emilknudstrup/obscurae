.. _install:

Installation
============

.. code-block:: bash
        
    cd /path/to
    git clone https://github.com/emilknudstrup/obscurae.git
    cd /path/to/fetch_n_prep
    python -m pip install .

or for an editable version

.. code-block:: bash

    python -m pip install -e .


Testing
============

You can test if the installation was successful by running the following command

.. code-block:: bash

    python -c 'import obscurae; obscurae.test()'