#include <stdio.h>
#include <string.h>
#include "animation.h"
#include "interpolation.h"

int main(int argc, char *argv[]) {
    if (argc < 2) {
        printf("Veuillez spécifier un argument : interpolation ou animation\n");
        return 1;
    }

    if (strcmp(argv[1], "interpolation") == 0) {
        run_interpolation();
    } else if (strcmp(argv[1], "animation") == 0) {
        run_animation();
    } else {
        printf("Argument inconnu : %s\n", argv[1]);
        printf("Utilisez 'interpolation' ou 'animation'\n");
        return 1;
    }

    return 0;
}
