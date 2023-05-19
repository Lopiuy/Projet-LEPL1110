<h1 align="center">
  <br>
    Projet Elasticité Linéaire
  <br>
</h1>

<h4 align="center">LEPL1110 - Eléments finis</h4>

<p align="center">
  <a href="#definition-du-problème">Définition du problème</a> •
  <a href="#usage">Usage</a> •
  <a href="#exécution">Exécution</a> •
  <a href="#auteurs">Auteurs</a> •
  <a href="#ressources">Ressources</a>
</p>

## Définition du problème

Le code proposé dans ce répertoire permet de résoudre différents problèmes d'élasticité linéaire en utilisant la méthode des éléments finis.
- On peut choisir le type de problèmes :
  - problème en déformation plane
  - problème en tension plane
  - problème axisymétrique

- On peut aussi choisir le type de résolution :
  - résolution par matrice pleine
  - résolution par matrice bande

- On peut choisir le type de renumérotation des noeuds :
  - renumérotation selon leurs coordonnées en X
  - renumérotation selon leurs coordonnées en Y
  - pas de renumérotation

- On peut choisir le type de maillage :
  - maillage en triangle
  - maillage en quadrilatère

Nous avons également créer 3 problèmes d'exemples :
- problème du traîneau, qui représente le problème explicité dans le rapport
- problème de la maison, qui représente une maison avec une force sur le toit
- problème de la poutre
- problème simple (point 3.0.1 du document note-elasticite, situé sur le site du cours)

## Usage

Les choix proposés ci-dessus peuvent être modifiés dans le fichier <em>src/ProjectPreProcessor/src/main.c</em> dans la fonction suivante:
```c 
femElasticityCreate(femGeo* theGeometry, double E, double nu, double rho, double g, femElasticCase iCase, femSolverType iSolver, femRenumType iRenum)
```
- le type de problème peut être modifié en modifiant la variable `iCase` dans la fonction <em>main</em>
  - `PLANAR_STRAIN` : deformation plane
  - `PLANAR_STRESS` : tension plane
  - `AXISYMMETRIC` : axisymétrique

- le type de résolution peut être modifié en modifiant la variable `iSolver` dans la fonction <em>main</em>
  - `FEM_FULL` : matrice pleine
  - `FEM_BAND` : matrice bande

- le type de renumérotation peut être modifié en modifiant la variable `iRenum` dans la fonction <em>main</em>
  - `FEM_NO` : pas de renumérotation
  - `FEM_XNUM` : renumérotation selon les coordonnées en X
  - `FEM_YNUM` : renumérotation selon les coordonnées en Y

Pour modifier le type de maillage, il faut modifier la variable 'theGeometry->elementType' dans le fichier <em>src/ProjectPreProcessor/src/main.c</em> dans la fonction <em>main</em>
  - `FEM_TRIANGLE` : maillage en triangle
  - `FEM_QUAD` : maillage en quadrilatère

Dans la fonction nous retrouvons également les constantes physiques suivantes:
- `E` : module d'Young
- `nu` : coefficient de Poisson
- `rho` : masse volumique
- `g` : gravité

Le choix de problèmes d'exemple est fait lors de l'exécution comme expliqué dans la section [Execution](#execution).


## Exécution

Les différentes commandes proposées ci-dessous sont à compiler sous Linux.
Les scripts proposés sont à éxécuter depuis le répertoire <em>src</em> du projet.

```bash
# Pour acceder au repertoire src
$ cd src
# Pour construire le projet sur votre machine
$ ./run.sh build
```
```bash
# Pour lancer la résolution de différents problèmes
$ ./run.sh [problème] [forme]
# - [problème] : maison, poutre, basic
#    * luge : résolution du problème du traîneau
#    * maison : résolution du problème de la maison
#    * poutre : résolution du problème de la poutre
#    * basic : résolution du problème simple
# - [forme] : form
#    * form : affiche la forme du problème avant résolution
```
>**Note** 
> - Le script <em>run.sh</em> peut aussi être executer sans argument pour résoudre notre problème 'luge' par défaut.
> - Les différents problèmes qu'on peut mettre en argument de *[problème]* peuvent être modifiés comme explicité dans la
> section [Usage](#usage).

## Auteurs

Ce répertoire a été entièrement implémenté par Ygor Lausberg et Antoine Lemaire dans le cadre du cours de *Element Finis* dispensé par Pr. Vincent Legat.

## Ressources

* [Eléments_finis - LEPL1110](https://perso.uclouvain.be/vincent.legat/zouLab/epl1110.php) - Site du cours
---

> GitHub [Projet](https://github.com/Lopiuy/Projet-LEPL1110) &nbsp;&middot;&nbsp;
> [Ygor Lausberg](mailto:ygor.lausberg@student.uclouvain.be) &nbsp;&middot;&nbsp;
> GitHub [@Lopiuy](https://github.com/Lopiuy) &nbsp;&middot;&nbsp;
> [Antoine Lemaire](mailto:antoine.g.lemaire@student.uclouvain.be) &nbsp;&middot;&nbsp;
> GitHub [@AntoineLem](https://github.com/AntoineLem)