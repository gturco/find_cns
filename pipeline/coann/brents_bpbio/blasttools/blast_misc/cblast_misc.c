#include <stdio.h>
#include <string.h>
#define MAXLINE 8192
#define DEBUG 1


int main(int argc, char *argv[]) {
    return add_locs(argv[1], argv[2]);
}

int add_locs(char *in_name, char *out_name){
    float pct_id, bitscore;
    char *tmp_name;

    if(out_name == NULL){
        char tmp[60];
        tmp_name = tmpnam(tmp);
    }
    else {
        tmp_name = out_name;
    }
    if(DEBUG){
        printf("name:%s\n", tmp_name);
    }

    float evalue;
    FILE *IN, *OUT;
    
    char query[128];
    char subject[128];
    unsigned int hit_length, nmiss, ngaps, qstart, qstop, sstart, sstop;
    unsigned int gchr, gstart, gstop;
    const char *in_format  = "%s\t%s\t%f\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%e\t%f\n";
    const char *out_format = "%s\t%s\t%2.2f\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%.2e\t%.1f\n";

    IN = fopen(in_name, "r");
    OUT = fopen(tmp_name, "w");
    while((fscanf(IN, in_format, query, subject, &pct_id, &hit_length, &nmiss
                , &ngaps , &qstart, &qstop, &sstart, &sstop, &evalue, &bitscore))
                != EOF){

        sscanf(query, "%d||%d||%d||", &gchr, &gstart, &gstop);
        --gstart; // starts at 1, not 0. 1 + 1 => 1.
        qstart += gstart;
        qstop  += gstart;

        sscanf(subject, "%d||%d||%d||", &gchr, &gstart, &gstop);
        --gstart; // starts at 1, not 0. 1 + 1 => 1.
        sstart += gstart;
        sstop  += gstart;

        fprintf(OUT, out_format, query, subject, pct_id, hit_length, nmiss, ngaps, qstart, qstop, sstart, sstop, evalue, bitscore);

    }
    fclose(IN);
    fclose(OUT);
    if(out_name == NULL){
        printf("OVERWRITING EXISTING BLAST FILE: %s WITH UPDATED LOCS.\n"
                , in_name);
        printf("code:%d\n", rename(tmp_name, in_name));
    }
    return 0;

}

//    const char *qq = "10||000090267||000101787||EVM PREDICTION SUPERCONTIG_10.4||-1||CDS||94616543";
// sscanf("10||000090267||000101787||EVM PREDICTION SUPERCONTIG_10.4||-1||CDS||94616543", "%d||%d||%d||%[^|]||%d||%3s||%d ", &gchr, &gstart, &gstop, gname, &gstrand, gtype, &gfeature_id);
