.. installation:

Installation
==============

Installation from source
--------------------------

.. code-block:: bash

    $ conda create --name molfrag python=3.11
    $ conda activate molfrag
    $ git clone https://github.com/ckz1/MolFragApp.git
    $ pip install -r requirements.txt

Test your installation
------------------------

- Stop firewall (Linux) :code:`systemctl stop firewalld`
- Run directly :code:`streamlit run MolFragApp.py` or run in background :code:`nohup streamlit run MolFragApp.py > MolFragApp.log 2>&1 &`