import time
import argparse
import sqlite3
import re
import multiprocessing
from collections import Counter

TIMEOUT_LOAD_SERVER = 1800

LEVEL_CONTENT = {
    "NOG": [
        "arNOG",
        "bactNOG",
        "euNOG",
        "thaNOG",
        "eurNOG",
        "creNOG",
        "synNOG",
        "spiNOG",
        "firmNOG",
        "fusoNOG",
        "aquNOG",
        "aciNOG",
        "therNOG",
        "tenNOG",
        "proNOG",
        "defNOG",
        "plaNOG",
        "actNOG",
        "chloNOG",
        "cyaNOG",
        "deiNOG",
        "bctoNOG",
        "chlNOG",
        "chlaNOG",
        "verNOG",
        "kinNOG",
        "perNOG",
        "virNOG",
        "apiNOG",
        "opiNOG",
        "arcNOG",
        "metNOG",
        "methNOG",
        "thermNOG",
        "methaNOG",
        "halNOG",
        "theNOG",
        "negNOG",
        "cloNOG",
        "eryNOG",
        "bacNOG",
        "acidNOG",
        "delNOG",
        "gproNOG",
        "aproNOG",
        "bproNOG",
        "chlorNOG",
        "dehNOG",
        "cytNOG",
        "bacteNOG",
        "sphNOG",
        "flaNOG",
        "verrNOG",
        "strNOG",
        "chloroNOG",
        "acoNOG",
        "cocNOG",
        "meNOG",
        "fuNOG",
        "dproNOG",
        "eproNOG",
        "braNOG",
        "lilNOG",
        "haeNOG",
        "cryNOG",
        "biNOG",
        "basNOG",
        "ascNOG",
        "poaNOG",
        "nemNOG",
        "artNOG",
        "chorNOG",
        "agarNOG",
        "treNOG",
        "saccNOG",
        "euroNOG",
        "sordNOG",
        "dotNOG",
        "chrNOG",
        "inNOG",
        "veNOG",
        "agaNOG",
        "sacNOG",
        "debNOG",
        "eurotNOG",
        "onyNOG",
        "hypNOG",
        "magNOG",
        "sorNOG",
        "pleNOG",
        "rhaNOG",
        "lepNOG",
        "dipNOG",
        "hymNOG",
        "fiNOG",
        "aveNOG",
        "maNOG",
        "arthNOG",
        "necNOG",
        "chaNOG",
        "droNOG",
        "spriNOG",
        "carNOG",
        "prNOG",
        "roNOG",
        "homNOG",
    ],
    "arNOG": [
        "thaNOG",
        "eurNOG",
        "creNOG",
        "arcNOG",
        "metNOG",
        "methNOG",
        "thermNOG",
        "methaNOG",
        "halNOG",
        "theNOG",
    ],
    "bactNOG": [
        "synNOG",
        "spiNOG",
        "firmNOG",
        "fusoNOG",
        "aquNOG",
        "aciNOG",
        "therNOG",
        "tenNOG",
        "proNOG",
        "defNOG",
        "plaNOG",
        "actNOG",
        "chloNOG",
        "cyaNOG",
        "deiNOG",
        "bctoNOG",
        "chlNOG",
        "chlaNOG",
        "verNOG",
        "negNOG",
        "cloNOG",
        "eryNOG",
        "bacNOG",
        "acidNOG",
        "delNOG",
        "gproNOG",
        "aproNOG",
        "bproNOG",
        "chlorNOG",
        "dehNOG",
        "cytNOG",
        "bacteNOG",
        "sphNOG",
        "flaNOG",
        "verrNOG",
        "dproNOG",
        "eproNOG",
    ],
    "euNOG": [
        "kinNOG",
        "perNOG",
        "virNOG",
        "apiNOG",
        "opiNOG",
        "strNOG",
        "chloroNOG",
        "acoNOG",
        "cocNOG",
        "meNOG",
        "fuNOG",
        "braNOG",
        "lilNOG",
        "haeNOG",
        "cryNOG",
        "biNOG",
        "basNOG",
        "ascNOG",
        "poaNOG",
        "nemNOG",
        "artNOG",
        "chorNOG",
        "agarNOG",
        "treNOG",
        "saccNOG",
        "euroNOG",
        "sordNOG",
        "dotNOG",
        "chrNOG",
        "inNOG",
        "veNOG",
        "agaNOG",
        "sacNOG",
        "debNOG",
        "eurotNOG",
        "onyNOG",
        "hypNOG",
        "magNOG",
        "sorNOG",
        "pleNOG",
        "rhaNOG",
        "lepNOG",
        "dipNOG",
        "hymNOG",
        "fiNOG",
        "aveNOG",
        "maNOG",
        "arthNOG",
        "necNOG",
        "chaNOG",
        "droNOG",
        "spriNOG",
        "carNOG",
        "prNOG",
        "roNOG",
        "homNOG",
    ],
    "eurNOG": [
        "arcNOG",
        "metNOG",
        "methNOG",
        "thermNOG",
        "methaNOG",
        "halNOG",
        "theNOG",
    ],
    "firmNOG": ["negNOG", "cloNOG", "eryNOG", "bacNOG"],
    "aciNOG": ["acidNOG"],
    "proNOG": ["delNOG", "gproNOG", "aproNOG", "bproNOG", "dproNOG", "eproNOG"],
    "chloNOG": ["chlorNOG", "dehNOG"],
    "bctoNOG": ["cytNOG", "bacteNOG", "sphNOG", "flaNOG"],
    "verNOG": ["verrNOG"],
    "virNOG": ["strNOG", "chloroNOG", "braNOG", "lilNOG", "poaNOG"],
    "apiNOG": ["acoNOG", "cocNOG", "haeNOG", "cryNOG"],
    "opiNOG": [
        "meNOG",
        "fuNOG",
        "biNOG",
        "basNOG",
        "ascNOG",
        "nemNOG",
        "artNOG",
        "chorNOG",
        "agarNOG",
        "treNOG",
        "saccNOG",
        "euroNOG",
        "sordNOG",
        "dotNOG",
        "chrNOG",
        "inNOG",
        "veNOG",
        "agaNOG",
        "sacNOG",
        "debNOG",
        "eurotNOG",
        "onyNOG",
        "hypNOG",
        "magNOG",
        "sorNOG",
        "pleNOG",
        "rhaNOG",
        "lepNOG",
        "dipNOG",
        "hymNOG",
        "fiNOG",
        "aveNOG",
        "maNOG",
        "arthNOG",
        "necNOG",
        "chaNOG",
        "droNOG",
        "spriNOG",
        "carNOG",
        "prNOG",
        "roNOG",
        "homNOG",
    ],
    "delNOG": ["dproNOG", "eproNOG"],
    "strNOG": ["braNOG", "lilNOG", "poaNOG"],
    "acoNOG": ["haeNOG"],
    "cocNOG": ["cryNOG"],
    "meNOG": [
        "biNOG",
        "nemNOG",
        "artNOG",
        "chorNOG",
        "chrNOG",
        "inNOG",
        "veNOG",
        "rhaNOG",
        "lepNOG",
        "dipNOG",
        "hymNOG",
        "fiNOG",
        "aveNOG",
        "maNOG",
        "droNOG",
        "spriNOG",
        "carNOG",
        "prNOG",
        "roNOG",
        "homNOG",
    ],
    "fuNOG": [
        "basNOG",
        "ascNOG",
        "agarNOG",
        "treNOG",
        "saccNOG",
        "euroNOG",
        "sordNOG",
        "dotNOG",
        "agaNOG",
        "sacNOG",
        "debNOG",
        "eurotNOG",
        "onyNOG",
        "hypNOG",
        "magNOG",
        "sorNOG",
        "pleNOG",
        "arthNOG",
        "necNOG",
        "chaNOG",
    ],
    "lilNOG": ["poaNOG"],
    "biNOG": [
        "nemNOG",
        "artNOG",
        "chorNOG",
        "chrNOG",
        "inNOG",
        "veNOG",
        "rhaNOG",
        "lepNOG",
        "dipNOG",
        "hymNOG",
        "fiNOG",
        "aveNOG",
        "maNOG",
        "droNOG",
        "spriNOG",
        "carNOG",
        "prNOG",
        "roNOG",
        "homNOG",
    ],
    "basNOG": ["agarNOG", "treNOG", "agaNOG"],
    "ascNOG": [
        "saccNOG",
        "euroNOG",
        "sordNOG",
        "dotNOG",
        "sacNOG",
        "debNOG",
        "eurotNOG",
        "onyNOG",
        "hypNOG",
        "magNOG",
        "sorNOG",
        "pleNOG",
        "arthNOG",
        "necNOG",
        "chaNOG",
    ],
    "nemNOG": ["chrNOG", "rhaNOG"],
    "artNOG": ["inNOG", "lepNOG", "dipNOG", "hymNOG", "droNOG"],
    "chorNOG": [
        "veNOG",
        "fiNOG",
        "aveNOG",
        "maNOG",
        "spriNOG",
        "carNOG",
        "prNOG",
        "roNOG",
        "homNOG",
    ],
    "agarNOG": ["agaNOG"],
    "saccNOG": ["sacNOG", "debNOG"],
    "euroNOG": ["eurotNOG", "onyNOG", "arthNOG"],
    "sordNOG": ["hypNOG", "magNOG", "sorNOG", "necNOG", "chaNOG"],
    "dotNOG": ["pleNOG"],
    "chrNOG": ["rhaNOG"],
    "inNOG": ["lepNOG", "dipNOG", "hymNOG", "droNOG"],
    "veNOG": [
        "fiNOG",
        "aveNOG",
        "maNOG",
        "spriNOG",
        "carNOG",
        "prNOG",
        "roNOG",
        "homNOG",
    ],
    "onyNOG": ["arthNOG"],
    "hypNOG": ["necNOG"],
    "sorNOG": ["chaNOG"],
    "dipNOG": ["droNOG"],
    "maNOG": ["spriNOG", "carNOG", "prNOG", "roNOG", "homNOG"],
    "spriNOG": ["prNOG", "roNOG", "homNOG"],
    "prNOG": ["homNOG"],
}


LEVEL_HIERARCHY = (
    "((thaNOG,(arcNOG,metNOG,methNOG,thermNOG,methaNOG,halNOG,theNOG)eurNOG,creNOG)arNOG,(synNOG,spiNOG,(negNOG,cloNOG,eryNOG,bacNOG)firmNOG,fusoNOG,aquNOG,(acidNOG)aciNOG,therNOG,tenNOG,((cytNOG,bacteNOG,sphNOG,flaNOG)bctoNOG,chlNOG)NoName,((dproNOG,eproNOG)delNOG,gproNOG,aproNOG,bproNOG)proNOG,defNOG,plaNOG,actNOG,(chlorNOG,dehNOG)chloNOG,cyaNOG,(chlaNOG,(verrNOG)verNOG)NoName,deiNOG)bactNOG,(kinNOG,perNOG,(((braNOG,(poaNOG)lilNOG)NoName)strNOG,chloroNOG)virNOG,((haeNOG)acoNOG,(cryNOG)cocNOG)apiNOG,(((((((lepNOG,(droNOG)dipNOG,hymNOG)NoName)inNOG)artNOG,(((fiNOG,(aveNOG,((((homNOG)prNOG,roNOG)spriNOG,carNOG)NoName)maNOG)NoName)NoName)veNOG)chorNOG)NoName,((rhaNOG)chrNOG)nemNOG)biNOG)meNOG,(((((agaNOG)agarNOG,treNOG)NoName)basNOG,(((((eurotNOG,(arthNOG)onyNOG)NoName)euroNOG,((magNOG,(chaNOG)sorNOG)NoName,(necNOG)hypNOG)sordNOG,(pleNOG)dotNOG)NoName,((sacNOG,debNOG)NoName)saccNOG)NoName)ascNOG)NoName)fuNOG)opiNOG)euNOG)NOG;"
)
EGGNOG_DATABASES = {
    k: 51700 + (i * 2)
    for i, k in enumerate(
        "NOG,aciNOG,acidNOG,acoNOG,actNOG,agaNOG,agarNOG,apiNOG,aproNOG,aquNOG,arNOG,arcNOG,artNOG,arthNOG,ascNOG,aveNOG,bacNOG,bactNOG,bacteNOG,basNOG,bctoNOG,biNOG,bproNOG,braNOG,carNOG,chaNOG,chlNOG,chlaNOG,chloNOG,chlorNOG,chloroNOG,chorNOG,chrNOG,cloNOG,cocNOG,creNOG,cryNOG,cyaNOG,cytNOG,debNOG,defNOG,dehNOG,deiNOG,delNOG,dipNOG,dotNOG,dproNOG,droNOG,eproNOG,eryNOG,euNOG,eurNOG,euroNOG,eurotNOG,fiNOG,firmNOG,flaNOG,fuNOG,fusoNOG,gproNOG,haeNOG,halNOG,homNOG,hymNOG,hypNOG,inNOG,kinNOG,lepNOG,lilNOG,maNOG,magNOG,meNOG,metNOG,methNOG,methaNOG,necNOG,negNOG,nemNOG,onyNOG,opiNOG,perNOG,plaNOG,pleNOG,poaNOG,prNOG,proNOG,rhaNOG,roNOG,sacNOG,saccNOG,sorNOG,sordNOG,sphNOG,spiNOG,spriNOG,strNOG,synNOG,tenNOG,thaNOG,theNOG,therNOG,thermNOG,treNOG,veNOG,verNOG,verrNOG,virNOG,viruses".split(
            ","
        )
    )
}
EGGNOG_DATABASES.update({"euk": 51400, "bact": 51500, "arch": 51600})
TAXID2LEVEL = {
    6656: "artNOG",
    1: "NOG",
    2: "bactNOG",
    5125: "hypNOG",
    5139: "sorNOG",
    5653: "kinNOG",
    110618: "necNOG",
    7711: "chorNOG",
    508458: "synNOG",
    639021: "magNOG",
    7214: "droNOG",
    28211: "aproNOG",
    28216: "bproNOG",
    28221: "dproNOG",
    7742: "veNOG",
    1090: "chlNOG",
    8782: "aveNOG",
    200783: "aquNOG",
    34384: "arthNOG",
    5204: "basNOG",
    147541: "dotNOG",
    6231: "nemNOG",
    147545: "euroNOG",
    200795: "chloNOG",
    6236: "rhaNOG",
    1117: "cyaNOG",
    147550: "sordNOG",
    92860: "pleNOG",
    909932: "negNOG",
    2157: "arNOG",
    5234: "treNOG",
    3699: "braNOG",
    183925: "metNOG",
    183939: "methNOG",
    204428: "chlaNOG",
    4751: "fuNOG",
    204432: "acidNOG",
    183963: "halNOG",
    183967: "theNOG",
    183968: "thermNOG",
    5794: "apiNOG",
    5796: "cocNOG",
    35493: "strNOG",
    4776: "perNOG",
    183980: "arcNOG",
    4893: "sacNOG",
    5819: "haeNOG",
    526524: "eryNOG",
    544448: "tenNOG",
    2759: "euNOG",
    1224: "proNOG",
    1236: "gproNOG",
    200918: "therNOG",
    1239: "firmNOG",
    28889: "creNOG",
    7898: "fiNOG",
    40674: "maNOG",
    9443: "prNOG",
    203494: "verrNOG",
    7399: "hymNOG",
    301297: "dehNOG",
    9989: "roNOG",
    35082: "cryNOG",
    1297: "deiNOG",
    33554: "carNOG",
    422676: "acoNOG",
    4890: "ascNOG",
    4891: "saccNOG",
    28890: "eurNOG",
    5338: "agaNOG",
    314146: "spriNOG",
    766764: "debNOG",
    119089: "chrNOG",
    32061: "chlorNOG",
    32066: "fusoNOG",
    200930: "defNOG",
    4447: "lilNOG",
    29547: "eproNOG",
    57723: "aciNOG",
    50557: "inNOG",
    651137: "thaNOG",
    33154: "opiNOG",
    9604: "homNOG",
    35718: "chaNOG",
    33090: "virNOG",
    33183: "onyNOG",
    203682: "plaNOG",
    38820: "poaNOG",
    203691: "spiNOG",
    68525: "delNOG",
    7088: "lepNOG",
    186801: "cloNOG",
    5042: "eurotNOG",
    91061: "bacNOG",
    33208: "meNOG",
    33213: "biNOG",
    200643: "bacteNOG",
    976: "bctoNOG",
    201174: "actNOG",
    74201: "verNOG",
    3041: "chloroNOG",
    155619: "agarNOG",
    7147: "dipNOG",
    117743: "flaNOG",
    117747: "sphNOG",
    224756: "methaNOG",
    768503: "cytNOG",
}
TAXONOMIC_RESOLUTION = [
    "apiNOG",
    "virNOG",
    "nemNOG",
    "artNOG",
    "maNOG",
    "fiNOG",
    "aveNOG",
    "meNOG",
    "fuNOG",
    "opiNOG",
    "euNOG",
    "arNOG",
    "bactNOG",
    "NOG",
]


def annotate_hits_file(seed_orthologs_file, annot_file, hmm_hits_file):
    annot_header = (
        "#query_name",
        "seed_eggNOG_ortholog",
        "seed_ortholog_evalue",
        "seed_ortholog_score",
        "predicted_gene_name",
        "GO_terms",
        "KEGG_KOs",
        "BiGG_reactions",
        "Annotation_tax_scope",
        "OGs",
        "bestOG|evalue|score",
        "COG cat",
        "eggNOG annot",
    )
    # logging.info('using annotate_hits_file')
    start_time = time.time()
    seq2bestOG = {}
    # if pexists(hmm_hits_file):
    seq2bestOG = get_seq_hmm_matches(hmm_hits_file)

    seq2annotOG = get_ogs_annotations(set([v[0] for v in seq2bestOG.values()]))

    print("Functional annotation of refined hits starts now", "green")
    OUT = open(annot_file, "w")
    # if args.report_orthologs:
    ORTHOLOGS = open(annot_file + ".orthologs", "w")

    # if not args.no_file_comments:
    # OUT.write('# emapper version:', get_version(), 'emapper DB:', get_db_version(), file=OUT)
    # print('# command: ./emapper.py ', ' '.join(sys.argv[1:]), file=OUT)
    # print('# time: ' + time.ctime(), file=OUT)
    OUT.write("\t".join(annot_header) + "\n")

    qn = 0
    # pool = multiprocessing.Pool(args.cpu)
    pool = multiprocessing.Pool(2)
    for result in pool.imap(annotate_hit_line, iter_hit_lines(seed_orthologs_file)):
        qn += 1
        if qn and (qn % 500 == 0):
            total_time = time.time() - start_time
            print(
                qn,
                total_time,
                "%0.2f q/s (refinement)" % ((float(qn) / total_time)),
                file=sys.stderr,
            )
            sys.stderr.flush()

        if result:
            (
                query_name,
                best_hit_name,
                best_hit_evalue,
                best_hit_score,
                best_name,
                gos,
                kegg,
                bigg,
                annot_level_max,
                match_nogs,
                orthologs,
            ) = result

            if query_name in seq2bestOG:
                (
                    hitname,
                    evalue,
                    score,
                    qlength,
                    hmmfrom,
                    hmmto,
                    seqfrom,
                    seqto,
                    q_coverage,
                ) = seq2bestOG[query_name]
                bestOG = "%s|%s|%s" % (hitname, evalue, score)
                og_cat, og_desc = seq2annotOG.get(hitname, ["", ""])
            else:
                bestOG = "NA|NA|NA"
                og_cat, og_desc = get_best_og_description(match_nogs)

            # if args.report_orthologs:
            ORTHOLOGS.write("\t".join(map(str, (query_name, ",".join(orthologs)))))

            OUT.write(
                "\t".join(
                    map(
                        str,
                        (
                            query_name,
                            best_hit_name,
                            best_hit_evalue,
                            best_hit_score,
                            best_name,
                            ",".join(sorted(gos)),
                            ",".join(sorted(kegg)),
                            ",".join(sorted(bigg)),
                            annot_level_max,
                            ",".join(match_nogs),
                            bestOG,
                            og_cat.replace("\n", ""),
                            og_desc.replace("\n", " "),
                        ),
                    )
                )
                + "\n"
            )

        # OUT.flush()
    pool.terminate()

    elapsed_time = time.time() - start_time
    # if not args.no_file_comments:
    # OUT.write('# %d queries scanned' % (qn))
    # OUT.write('# Total time (seconds):')
    # OUT.write('# Rate:', "%0.2f q/s" % ((float(qn) / elapsed_time)))
    OUT.close()

    # if args.report_orthologs:
    ORTHOLOGS.close()

    print(
        " Processed queries:%s total_time:%s rate:%s"
        % (qn, elapsed_time, "%0.2f q/s" % ((float(qn) / elapsed_time))),
        "lblue",
    )


def annotate_hit_line(arguments):
    # logging.info('using annotate_hit_line')
    connect()
    line = arguments

    if not line.strip() or line.startswith("#"):
        return None
    r = list(map(str.strip, line.split("\t")))

    query_name = r[0]
    best_hit_name = r[1]
    if best_hit_name == "-" or best_hit_name == "ERROR":
        return None

    best_hit_evalue = float(r[2])
    best_hit_score = float(r[3])
    if best_hit_score < 60 or best_hit_evalue > 0.001:
        return None

    match_nogs = get_member_ogs(best_hit_name)
    if not match_nogs:
        return None

    match_levels = set([nog.split("@")[1] for nog in match_nogs])
    # if args.tax_scope == "auto":
    for level in TAXONOMIC_RESOLUTION:
        if level in match_levels:
            annot_levels = set(LEVEL_CONTENT.get(level, [level]))
            annot_levels.add(level)
            annot_level_max = "%s[%d]" % (level, len(annot_levels))
            break
    # else:
    #     annot_levels = set(LEVEL_CONTENT.get(args.tax_scope, [args.tax_scope]))
    #     annot_levels.add(args.tax_scope)
    #     annot_level_max = "%s[%d]" %(args.tax_scope, len(annot_levels))

    all_orthologies = get_member_orthologs(best_hit_name, target_levels=annot_levels)
    ## just picking the target_orhtolog to be all. can change later
    # orthologs = sorted(all_orthologies[args.target_orthologs])
    orthologs = sorted(all_orthologies["all"])
    # if args.excluded_taxa:
    #     orthologs = [o for o in orthologs if not o.startswith("%s." %args.excluded_taxa)]

    if orthologs:
        pname, gos, kegg, bigg = summarize_annotations(
            orthologs, target_go_ev="non-electronic", excluded_go_ev=set(["ND", "IEA"])
        )

        best_name = ""
        if pname:
            name_candidate, freq = pname.most_common(1)[0]
            if freq >= 2:
                best_name = name_candidate
    else:
        pname = []
        best_name = ""
        gos = set()
        kegg = set()
        bigg = set()

    return (
        query_name,
        best_hit_name,
        best_hit_evalue,
        best_hit_score,
        best_name,
        gos,
        kegg,
        bigg,
        annot_level_max,
        match_nogs,
        orthologs,
    )


def summarize_annotations(seq_names, target_go_ev, excluded_go_ev):
    print("using summarize_annotations")
    in_clause = ",".join(['"%s"' % n for n in seq_names])
    cmd = (
        """SELECT eggnog.name, seq.pname, gene_ontology.gos, kegg.ko, bigg.reaction
        FROM eggnog
        LEFT JOIN seq on seq.name = eggnog.name
        LEFT JOIN gene_ontology on gene_ontology.name = eggnog.name
        LEFT JOIN kegg on kegg.name = eggnog.name
        LEFT JOIN bigg on bigg.name = eggnog.name
        WHERE eggnog.name in (%s)
        """
        % in_clause
    )
    all_gos = set()
    all_kegg = set()
    all_bigg = set()
    all_pnames = Counter()
    db.execute(cmd)
    for name, pname, gos, kegg, bigg in db.fetchall():
        if gos:
            all_gos.update(parse_gos(gos, target_go_ev, excluded_go_ev))
        if kegg:
            all_kegg.update([str(x).strip() for x in kegg.split(",")])
        if bigg:
            all_bigg.update([str(x).strip() for x in bigg.split(",")])
        if pname:
            all_pnames.update([pname.strip()])

    all_gos.discard("")
    all_kegg.discard("")
    all_bigg.discard("")
    del all_pnames[""]

    return all_pnames, all_gos, all_kegg, all_bigg


def parse_gos(gos, target_go_ev, excluded_go_ev):
    selected_gos = set()
    for g in gos.strip().split(","):
        if not g:
            continue
        gocat, gid, gevidence = list(map(str, g.strip().split("|")))
        if not target_go_ev or gevidence in target_go_ev:
            if not excluded_go_ev or gevidence not in excluded_go_ev:
                selected_gos.add(gid)
    return selected_gos


def get_member_ogs(name):
    print("using get_member_ogs")
    cmd = 'SELECT groups FROM eggnog WHERE name == "%s";' % (name)
    db.execute(cmd)
    match = db.fetchone()
    ogs = None
    if match:
        ogs = [str(x).strip() for x in match[0].split(",")]
    return ogs


def connect():
    print("using connect")
    global conn, db
    conn = sqlite3.connect(
        "/Users/zavo603/anaconda3/envs/python2/lib/python2.7/site-packages/data/eggnog.db"
    )
    db = conn.cursor()


def get_seq_hmm_matches(hits_file):
    # loggin.info('using get_seq_hmm_matches')
    connect()
    print("Reading HMM matches", "green")
    seq2oginfo = {}
    start_time = time.time()
    hitnames = set()
    # if pexists(hits_file):
    for line in open(hits_file):
        if not line.strip() or line.startswith("#"):
            continue
        print(line.strip().split())
        (
            query,
            hit,
            evalue,
            sum_score,
            query_length,
            hmmfrom,
            hmmto,
            seqfrom,
            seqto,
            q_coverage,
        ) = list(map(str.strip, line.split()))

        if query not in seq2oginfo and hit not in ["ERROR", "-"]:
            hitname = cleanup_og_name(hit)
            seq2oginfo[query] = [
                hitname,
                evalue,
                sum_score,
                query_length,
                hmmfrom,
                hmmto,
                seqfrom,
                seqto,
                q_coverage,
            ]
    return seq2oginfo


def get_member_orthologs(member, target_taxa=None, target_levels=None):
    print("using get_member_orthologs")
    query_taxa = member.split(".", 1)[0]
    target_members = set([member])
    cmd = 'SELECT orthoindex FROM orthologs WHERE name = "%s"' % member.strip()
    db.execute(cmd)
    event_indexes = str(db.fetchone()[0])
    cmd2 = "SELECT level, side1, side2 FROM event WHERE i IN (%s)" % event_indexes
    if target_levels:
        cmd2 += " AND level IN (%s)" % (",".join(['"%s"' % x for x in target_levels]))
    db.execute(cmd2)
    orthology = {}
    for level, _side1, _side2 in db.fetchall():
        side1 = [m.split(".", 1) for m in _side1.split(",")]
        side2 = [m.split(".", 1) for m in _side2.split(",")]

        by_sp1, by_sp2 = {}, {}
        for _sp, _side in [(by_sp1, side1), (by_sp2, side2)]:
            for t, s in _side:
                if not target_taxa or t in target_taxa or t == query_taxa:
                    mid = "%s.%s" % (t, s)
                    _sp.setdefault(t, set()).add(mid)

        # merge by side1 coorthologs
        targets = target_taxa or list(by_sp2.keys())
        for sp1, co1 in by_sp1.items():
            if target_members & co1:
                key1 = (sp1, tuple(sorted((co1))))
                for sp2 in targets:
                    if sp2 not in by_sp2:
                        continue
                    co2 = by_sp2[sp2]
                    key2 = (sp2, tuple(sorted(co2)))
                    orthology.setdefault(key1, set()).add(key2)

        # merge by side2 coorthologs
        targets = target_taxa or list(by_sp1.keys())
        for sp1, co1 in by_sp2.items():
            if target_members & co1:
                key1 = (sp1, tuple(sorted((co1))))
                for sp2 in targets:
                    if sp2 not in by_sp1:
                        continue
                    co2 = by_sp1[sp2]
                    key2 = (sp2, tuple(sorted(co2)))
                    orthology.setdefault(key1, set()).add(key2)

    all_orthologs = {
        "one2one": set(),
        "one2many": set(),
        "many2many": set(),
        "many2one": set(),
        "all": set(),
    }

    for k, v in orthology.items():
        if len(k[1]) == 1:
            otype_prefix = "one2"
        else:
            otype_prefix = "many2"
        all_orthologs["all"].update(k[1])
        for t2, co2 in v:
            if len(co2) == 1:
                otype = otype_prefix + "one"
            else:
                otype = otype_prefix + "many"
            all_orthologs[otype].update(k[1])
            all_orthologs[otype].update(co2)
            all_orthologs["all"].update(co2)

    return all_orthologs


def get_ogs_annotations(ognames):
    print("using get_ogs_annotations")
    # og VARCHAR(16) PRIMARY KEY,
    # level VARCHAR(16),
    # nm INTEGER,
    # description TEXT,
    # COG_categories VARCHAR(8),
    # GO_freq TEXT,
    # KEGG_freq TEXT,
    # SMART_freq TEXT,
    # proteins TEXT);
    query = ",".join(['"%s"' % x for x in ognames])
    cmd = (
        "SELECT og.og, description, COG_categories FROM og WHERE og.og IN (%s)" % query
    )
    og2desc = {}
    if db.execute(cmd):
        cog_cat_cleaner = re.compile("[\[\]u'\"]+")
        for og, desc, cat in db.fetchall():
            cat = re.sub(cog_cat_cleaner, "", cat)
            og2desc[og] = [cat, desc]
    return og2desc


def iter_hit_lines(filename):
    # logging.info('using iter_hit_lines')
    for line in open(filename):
        if line.startswith("#") or not line.strip():
            continue
        yield line


def cleanup_og_name(name):
    # names in the hmm databases are sometiemes not clean eggnog OG names
    # logging.info('using cleanup_og_name')
    m = re.search("\w+\.((ENOG41|COG|KOG|arCOG)\w+)\.", name)
    if m:
        name = m.groups()[0]
    name = re.sub("^ENOG41", "", name)
    return name


def main(search_output, annot_file, hits_file):
    annotate_hits_file(search_output, annot_file, hits_file)


if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("--search_output", type=str, help="output file of seq_search.py")
    p.add_argument("anot_file", type=str, help="annotation ouptut file name")
    p.add_argument("--hits_file", type=str, help="hmm hits file returned from hmmscan")

    args = p.parse_args()
    # logging.basicConfig(
    #     level=logging.INFO, datefmt="%Y-%m-%d %H:%M", format="[%(asctime)s] %(message)s"
    # )
    main(args.search_output, args.anot_file, args.hits_file)
