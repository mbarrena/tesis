pip install matplotlib
pip install seaborn
pip install --upgrade 'xlrd>=1.2.0'
pip install scikit-learn
pip install ruptures
#pip install chowtest   #Si falla esta lista ir a Entorno de ejecucion > reiniciar entorno de ejecucion y volver a correr esta celda.
pip3 uninstall statsmodels -y
pip install statsmodels --upgrade
pip install arch
rm -rf chowtest
git clone https://github.com/mbarrena/chowtest.git && cd chowtest && pip install .
git clone https://github.com/mbarrena/localprojections.git && cd localprojections && pip install .
pip install --upgrade nbformat
