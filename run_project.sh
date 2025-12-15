#!/bin/bash

# Vérification si le script est exécuté avec les droits sudo
if [ "$(id -u)" -ne 0 ]; then
    echo "Ce script doit être exécuté avec des privilèges sudo." 
    exit 1
fi

# Mettre à jour les informations sur les paquets et installer les bibliothèques nécessaires
echo "Mise à jour des paquets et installation des bibliothèques essentielles..."
sudo apt-get update -y
sudo apt-get install -y build-essential

# Vérification si l'installation de gcc et make a réussi
if ! command -v gcc &> /dev/null; then
    echo "Erreur : gcc n'est pas installé. Veuillez installer build-essential."
    exit 1
fi

# Vérification si le répertoire results existe, sinon le créer
if [ ! -d "results" ]; then
    echo "Création du répertoire results..."
    mkdir results
fi

# Compilation du programme C
echo "Compilation du programme C..."
gcc -O2 -o results/output src/C_code.c -lm
if [ $? -ne 0 ]; then
    echo "Erreur lors de la compilation du programme C."
    exit 1
fi

# Vérification si l'exécutable existe et est exécutable
if [ ! -x "results/output" ]; then
    echo "Erreur : l'exécutable 'results/output' n'existe pas ou n'est pas exécutable."
    exit 1
fi

# Exécution du programme C
echo "Exécution du programme C..."
./results/output
if [ $? -ne 0 ]; then
    echo "Erreur lors de l'exécution du programme C."
    exit 1
fi

# Vérification si MATLAB est installé
echo "Vérification de l'installation de MATLAB..."
if ! command -v /usr/local/bin/matlab-2021b &> /dev/null; then
    echo "Erreur : MATLAB n'est pas installé ou n'est pas dans le chemin spécifié."
    exit 1
fi

# Exécution du script MATLAB
echo "Exécution du script MATLAB..."
/usr/local/bin/matlab-2021b -batch "run('src/Matlab_code.m')"
if [ $? -ne 0 ]; then
    echo "Erreur lors de l'exécution du script MATLAB."
    exit 1
fi

echo "Script terminé avec succès."
