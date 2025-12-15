#!/bin/bash

# Définir un répertoire local pour les bibliothèques
LIB_DIR="$HOME/libs"

# Vérification de l'existence du répertoire 'libs'
if [ ! -d "$LIB_DIR" ]; then
    echo "Création du répertoire pour les bibliothèques : $LIB_DIR"
    mkdir -p "$LIB_DIR"
fi

# Télécharger et installer GCC et les bibliothèques essentielles dans le répertoire local
echo "Téléchargement des sources nécessaires (si elles ne sont pas déjà présentes)..."

# Vérification si gcc est déjà présent localement, sinon installation
if ! command -v gcc &> /dev/null; then
    echo "GCC n'est pas installé localement. Téléchargement et installation..."

    # Télécharger GCC (et les outils nécessaires) depuis une source externe (par exemple, utiliser le tarball précompilé)
    wget -q https://ftp.gnu.org/gnu/gcc/gcc-10.2.0/gcc-10.2.0.tar.gz -P "$LIB_DIR"
    cd "$LIB_DIR" || exit
    tar -xzf gcc-10.2.0.tar.gz

    # Compiling GCC (une version minimale)
    cd gcc-10.2.0 || exit
    ./configure --prefix="$LIB_DIR/gcc-install" --disable-multilib
    make -j$(nproc)
    make install
fi

# Assurez-vous que le chemin du compilateur gcc local est correctement ajouté au PATH
export PATH="$LIB_DIR/gcc-install/bin:$PATH"

# Télécharger et installer les bibliothèques nécessaires (libc, libm, etc.) dans le répertoire local
echo "Téléchargement et installation des bibliothèques de base (libm, etc.)..."

# Télécharger la bibliothèque mathématique et les en-têtes de développement standard
# Nous allons supposer que libm est déjà incluse dans GCC, mais dans un autre cas tu pourrais le faire à partir des tarballs

# Télécharger les sources de GCC si nécessaire
if [ ! -d "$LIB_DIR/gcc-install" ]; then
    echo "Les bibliothèques de GCC ne sont pas installées, installation en cours."
    cd "$LIB_DIR" || exit
    wget -q https://ftp.gnu.org/gnu/gcc/gcc-10.2.0/gcc-10.2.0.tar.gz
    tar -xzf gcc-10.2.0.tar.gz
    cd gcc-10.2.0 || exit
    ./configure --prefix="$LIB_DI_"
