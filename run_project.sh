#!/bin/bash

# Définir un répertoire local pour les bibliothèques
LIB_DIR="$HOME/libs"

# Vérification de l'existence du répertoire 'libs'
if [ ! -d "$LIB_DIR" ]; then
    echo "Création du répertoire pour les bibliothèques : $LIB_DIR"
    mkdir -p "$LIB_DIR"
fi

# Vérification de l'existence des bibliothèques nécessaires (GMP, MPFR, MPC) et installation locale
echo "Vérification des bibliothèques nécessaires (GMP, MPFR, MPC)..."

# Télécharger et installer GMP, MPFR et MPC si nécessaire

# Compilation de GCC : en tenant compte de l'option --with-gmp, --with-mpfr et --with-mpc
echo "Compilation de GCC (si nécessaire)..."
cd "$LIB_DIR" || exit
wget -q https://ftp.gnu.org/gnu/gcc/gcc-10.2.0/gcc-10.2.0.tar.gz
tar -xzf gcc-10.2.0.tar.gz
cd gcc-10.2.0 || exit

# Compiler GCC avec les bibliothèques locales
./configure --prefix="$LIB_DIR/gcc-install" --with-gmp="$LIB_DIR/gmp" --with-mpfr="$LIB_DIR/mpfr" --with-mpc="$LIB_DIR/mpc" --disable-multilib
make -j$(nproc)
make install

# Ajouter GCC local dans le PATH
export PATH="$LIB_DIR/gcc-install/bin:$PATH"

# Vérification si gcc est bien installé localement
if ! command -v gcc &> /dev/null; then
    echo "Erreur : gcc n'est pas installé localement."
    exit 1
fi

# Vérification si le répertoire results existe
if [ ! -d "results" ]; then
    echo "Création du répertoire 'results'..."
    mkdir results
fi

# Vérifier si les fichiers CSV existent déjà, sinon on les génère avec le programme C
if [ ! -f "~/CMT-Project/results/healthy_results.csv" ] || [ ! -f "~/CMT-Project/results/smoker_results.csv" ]; then
    echo "Les fichiers CSV sont manquants, exécution du programme C pour les générer..."
    
    # Compiler le programme C
    echo "Compilation du programme C..."
    gcc -O2 -o results/output ~/CMT-Project/src/C_code.c -lm
    if [ $? -ne 0 ]; then
        echo "Erreur lors de la compilation du programme C."
        exit 1
    fi

    # Exécuter le programme C pour générer les fichiers CSV
    echo "Exécution du programme C..."
    ./results/output
    if [ $? -ne 0 ]; then
        echo "Erreur lors de l'exécution du programme C."
        exit 1
    fi
else
    echo "Les fichiers CSV existent déjà, pas besoin de les régénérer."
fi

# Exécution du script MATLAB (assurer que le chemin vers le fichier Matlab est correct)
echo "Exécution du script MATLAB..."
/usr/local/bin/matlab-2021b -batch "run('$HOME/CMT-Project/src/Matlab_code.m')"
if [ $? -ne 0 ]; then
    echo "Erreur lors de l'exécution du script MATLAB."
    exit 1
fi

echo "Script terminé avec succès."
