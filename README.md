# Projet LINMA1100 : Analyse de la déformation d'un barrage hydraulique

Ce projet, réalisé par **Jean DE WALQUE** et **Gaston DE CUMONT**, vise à analyser la déformation d'une section horizontale d'un barrage hydraulique sous l'effet de la pression hydrostatique de l'eau.

## Structure du Projet

- **`build/`** : Contient les fichiers générés par CMake et les binaires.
- **`data/`** : Contient les fichiers de données nécessaires pour les calculs (ex. : maillages, paramètres d'élasticité).
- **`glfw/`** : Contient les fichiers liés à la bibliothèque GLFW.
- **`gmsh/`** : Contient le SDK de GMSH pour la génération de maillages.
- **`src/`** : Contient le code source principal du projet.

### Contenu du dossier `src/`

Le dossier `src/` contient les fichiers source nécessaires à l'exécution du projet. Voici une description détaillée des principaux fichiers :

- **`main.c`** : Le point d'entrée du programme. Ce fichier gère les arguments passés en ligne de commande (`interpolation` ou `animation`) et appelle les fonctions correspondantes.
- **`animation.h` et `animation.c`** : Contiennent les fonctions nécessaires pour gérer et afficher l'animation de la déformation du barrage.
- **`interpolation.h` et `interpolation.c`** : Implémentent les calculs d'interpolation, notamment l'interpolation d'Hermite utilisée pour ajuster les maillages.
- **`flag.h`** : Déclare la variable globale `flag`, utilisée pour différencier les modes d'exécution (`interpolation` ou `animation`). Car nous ne faisons pas d'interpollation d'hermite pour le maillage de notre barrage lors de l'animation.
- **`utils.h` et `utils.c`** : Fournissent des fonctions utilitaires, comme la génération de maillages (`geoMeshGenerate`) et des calculs géométriques (ex. : `geoSize`).
- **`fem.h` et `fem.c`** : Implémentent les fonctions liées aux éléments finis, comme l'assemblage des matrices et la résolution des systèmes linéaires.


## Prérequis

Avant de compiler et d'exécuter le projet, assurez-vous d'avoir installé les outils suivants :

- **CMake** : Pour la configuration du projet.
- **Compilateur C** : GCC ou équivalent.
- **GLFW** : Pour les animations graphiques.
- **GMSH** : Pour la génération de maillages.

## Compilation

1. Créez et entrez dans le répertoire `build` :
   ```bash
   mkdir build && cd build
   ```
2. Configurez le projet avec CMake :
   ```bash
   cmake ..
   ```
3. Compilez le programme :
   ```bash
   make
   ```

## Exécution

Une fois le programme compilé, vous pouvez l'exécuter avec les arguments suivants :

- **Pour analyser la déformation du barrage à une hauteur fixe (H = 50 m), en utilisant une interpollation d'hermite pour le maillage** :
  ```bash
  ./myFem interpolation
  ```

- **Pour visualiser une animation de la déformation du barrage avec une hauteur variant de 1 m à 50 m** :
  ```bash
  ./myFem animation
  ```

Si aucun argument ou un argument invalide est fourni, le programme affichera un message d'erreur expliquant les options disponibles.

## Environnement de Travail

Nous avons chacun travaillé dans un environnement Linux (Ubuntu) via **Windows Subsystem for Linux (WSL)**.
## Auteurs

- **Jean DE WALQUE (53202200)**
- **Gaston DE CUMONT (64062200)**
