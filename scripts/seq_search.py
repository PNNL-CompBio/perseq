import os
import gzip
import multiprocessing
import argparse
import time
import re
from collections import defaultdict
import uuid
import subprocess
from tempfile import NamedTemporaryFile


def dump_hmm_matches(hmmscan_output, hits_file):
    hits_header = ("#query_name", "hit", "evalue", "sum_score", "query_length",
                   "hmmfrom", "hmmto", "seqfrom", "seqto", "query_coverage")
    # Cache previous results if resuming is enabled
    VISITED = set()
    #with open(output_file,'w') as OUT:
    OUT = open(hits_file, 'w')
    OUT.write('# ' + '\t'.join(hits_header)+'\n')

    total_time = 0
    last_time = time.time()
    start_time = time.time()
    qn = 0 # in case nothing to loop bellow
    # parse_hmmscan
    for qn, (name, elapsed, hits, querylen, seq) in enumerate(parse_hmmscan(hmmscan_output)):

        if elapsed == -1:
            # error occurred
            OUT.write('\t'.join([name] + ['ERROR'] * (len(hits_header) - 1))+'\n')
            # print('\t'.join(
            #     [name] + ['ERROR'] * (len(hits_header) - 1)), file=OUT)
        elif not hits:
            OUT.write('\t'.join([name] + ['-'] * (len(hits_header) - 1))+'\n')
        else:
            for hitindex, (hid, heval, hscore, hmmfrom, hmmto, sqfrom, sqto, domscore) in enumerate(hits):
                hitname = hid
                # if idmap:
                #     hitname = idmap[hid][0]

                OUT.write('\t'.join(map(str, [name, hitname, heval, hscore,
                                                 int(querylen), int(hmmfrom),
                                                 int(hmmto), int(sqfrom),
                                                 int(sqto),
                                                 float(sqto - sqfrom) / querylen]))+'\n')
        # OUT.flush()

        # monitoring
        total_time += time.time() - last_time
        last_time = time.time()
        if qn and (qn % 25 == 0):
            print(qn + \
                1, total_time, "%0.2f q/s" % (float(qn + 1) / total_time))
            #sys.stderr.flush()

    # Writes final stats
    elapsed_time = time.time() - start_time
    #if not args.no_file_comments:
    # below go to OUT(in original)
    print('# %d queries scanned' % (qn + 1))
    print('# Total time (seconds):', elapsed_time)
    print('# Rate:', "%0.2f q/s" % ((float(qn + 1) / elapsed_time)))
    OUT.close()
    print(" Processed queries:%s total_time:%s rate:%s" %\
                   (qn+1, elapsed_time, "%0.2f q/s" % ((float(qn+1) / elapsed_time))), 'lblue')



# def iter_hits(source, translate, query_type, dbtype, scantype, host, port,
#               evalue_thr=None, score_thr=None, max_hits=None, return_seq=False,
#               skip=None, maxseqlen=None, fixed_Z=None, qcov_thr=None, cpus=1,
#               base_tempdir=None):
#     print('using iter_hits')
#     try:
#         max_hits = int(max_hits)
#     except Exception:
#         max_hits = None
#
#     # if scantype == 'mem' and query_type == "seq":
#     #     return iter_seq_hits(source, translate, host, port, dbtype=dbtype, evalue_thr=evalue_thr, score_thr=score_thr, max_hits=max_hits, skip=skip, maxseqlen=maxseqlen)
#     # elif scantype == 'mem' and query_type == "hmm" and dbtype == "seqdb":
#     #     return iter_hmm_hits(source, host, port, maxseqlen=maxseqlen)
#     return hmmscan(source, translate, host, evalue_thr=evalue_thr, score_thr=score_thr, max_hits=max_hits, cpus=cpus, maxseqlen=maxseqlen, base_tempdir=base_tempdir)
#     # else:

def parse_hmmscan(hmmscan_output,evalue_thr=None,score_thr=60,max_hits=None):
    # print('using hmmscan')
    # if not HMMSCAN:
    #     raise ValueError('hmmscan not found in path')
    #
    # tempdir = mkdtemp(prefix='emappertmp_hmmscan_', dir=base_tempdir)
    #
    # OUT = NamedTemporaryFile(dir=tempdir, mode='w+t')
    # if translate or maxseqlen:
    #     if translate:
    #         print('translating query input file')
    #     Q = NamedTemporaryFile(mode = 'w+t')
    #     for name, seq in seqio.iter_fasta_seqs(query_file, translate=translate):
    #         if maxseqlen is None or len(seq) <= maxseqlen:
    #             print(">%s\n%s" % (name, seq), file=Q)
    #     Q.flush()
    #     query_file = Q.name
    #
    # cmd = '%s --cpu %s -o /dev/null --domtblout %s %s %s' % (
    #     HMMSCAN, cpus, OUT.name, database_path, query_file)
    # print '#', cmd
    # print cmd
    #sts = subprocess.call(cmd, shell=True)
    byquery = defaultdict(list)

    last_query = None
    last_hitname = None
    hit_list = []
    hit_ids = set()
    last_query_len = None
    #if sts == 0:
    with open(hmmscan_output) as file:
        for line in file:
            # TBLOUT
            # ['#', '---', 'full', 'sequence', '----', '---', 'best', '1', 'domain', '----', '---', 'domain', 'number', 'estimation', '----']
            # ['#', 'target', 'name', 'accession', 'query', 'name', 'accession', 'E-value', 'score', 'bias', 'E-value', 'score', 'bias', 'exp', 'reg', 'clu', 'ov', 'env', 'dom', 'rep', 'inc', 'description', 'of', 'target']
            # ['#-------------------', '----------', '--------------------', '----------', '---------', '------', '-----', '---------', '------', '-----', '---', '---', '---', '---', '---', '---', '---', '---', '---------------------']
            # ['delNOG20504', '-', '553220', '-', '1.3e-116', '382.9', '6.2', '3.4e-116', '381.6', '6.2', '1.6', '1', '1', '0', '1', '1', '1', '1', '-']
            # fields = line.split() # output is not tab delimited! Should I trust this split?
            # hit, _, query, _ , evalue, score, bias, devalue, dscore, dbias = fields[0:10]

            # DOMTBLOUT
            #                                                                             --- full sequence --- -------------- this domain -------------   hmm coord   ali coord   env coord
            # target name        accession   tlen query name            accession   qlen   E-value  score  bias   #  of  c-Evalue  i-Evalue  score  bias  from    to  from    to  from    to  acc description of target
            # ------------------- ---------- -----  -------------------- -------
            # Pkinase              PF00069.22   264 1000565.METUNv1_02451 -
            # 858   4.5e-53  180.2   0.0   1   1   2.4e-56   6.6e-53  179.6
            # 0.0     1   253   580   830   580   838 0.89 Protein kinase
            # domain
            if line.startswith('#'):
                continue
            fields = line.split()

            (hitname, hacc, tlen, qname, qacc, qlen, evalue, score, bias, didx,
             dnum, c_evalue, i_evalue, d_score, d_bias, hmmfrom, hmmto, seqfrom,
             seqto, env_from, env_to, acc) = list(map(safe_cast, fields[:22]))

            if (last_query and qname != last_query):
                yield last_query, 0, hit_list, last_query_len, None
                hit_list = []
                hit_ids = set()
                last_query = qname
                last_query_len = None

            last_query = qname
            if last_query_len and last_query_len != qlen:
                raise ValueError(
                    "Inconsistent qlen when parsing hmmscan output")
            last_query_len = qlen

            if (evalue_thr is None or evalue <= evalue_thr) and \
               (score_thr is not None and score >= score_thr) and \
               (max_hits is None or last_hitname == hitname or len(hit_ids) < max_hits):

                hit_list.append([hitname, evalue, score, hmmfrom,
                                 hmmto, seqfrom, seqto, d_score])
                hit_ids.add(hitname)
                last_hitname = hitname

        if last_query:
            yield last_query, 0, hit_list, last_query_len, None

    # OUT.close()
    # if translate:
    #     Q.close()
    # shutil.rmtree(tempdir)
## gets used
def safe_cast(v):
    try:
        return float(v)
    except ValueError:
        return v.strip()

def refine_matches(fasta_file, refined_file, parsed_hmm_output,og2levels):
    refine_header = list(map(str.strip, '''#query_name, best_hit_eggNOG_ortholog,
                        best_hit_evalue, best_hit_score'''.split(',')))
    print('use refine_matches')
    print("Hit refinement starts now", 'green')
    start_time = time.time()
    #print(get_oglevels_file(),type(get_oglevels_file()))
    with gzip.open(og2levels,'rt') as file, open(refined_file,'w') as OUT:
        og2level = dict([tuple(map(str.strip, line.split('\t'))) for line in file])
    # og2level = dict([tuple(map(str.strip, line.split('\t')))
    #                  for line in gopen(get_oglevels_file())])
    #OUT = open(refine_file, "w")

    #if not args.no_file_comments:
    # below orignially written to out
        #OUT.write(get_call_info())
        OUT.write('\t'.join(refine_header))

        qn = 0 # in case no hits in loop below
        for qn, r in enumerate(process_nog_hits_file(parsed_hmm_output, fasta_file, og2level,
                                                     translate=False,
                                                     cpu=2,
                                                     excluded_taxa=None,
                                                     base_tempdir=None)):
            if qn and (qn % 25 == 0):
                total_time = time.time() - start_time
                print(qn + 1, total_time, "%0.2f q/s (refinement)" % ((float(qn + 1) / total_time)))
                #sys.stderr.flush()
            query_name = r[0]
            best_hit_name = r[1]
            if best_hit_name == '-' or best_hit_name == 'ERROR':
                continue
            best_hit_evalue = float(r[2])
            best_hit_score = float(r[3])
            OUT.write('\t'.join(map(str, (query_name, best_hit_name,
                                             best_hit_evalue, best_hit_score)))+'\n')
            #OUT.flush()

        elapsed_time = time.time() - start_time
        #if not args.no_file_comments:
            #below written to OUT
        print('# %d queries scanned' % (qn + 1))
        print('# Total time (seconds):', elapsed_time)
        print('# Rate:', "%0.2f q/s" % ((float(qn + 1) / elapsed_time)))
        #OUT.close()
        print(" Processed queries:%s total_time:%s rate:%s" %\
                       (qn+1, elapsed_time, "%0.2f q/s" % ((float(qn+1) / elapsed_time))), 'lblue')

def parse_seqs(translated_fasta):
    with open(translated_fasta) as file:
        n=0
        seq_dict={}
        for line in file:
            line = line.strip()
            if line.startswith('>'):
                header = line[1:].split()[0].strip()
            else:
                seq_dict[header]=line
    return seq_dict

def process_nog_hits_file(parsed_hmm_output, query_fasta, og2level, skip_queries=None,
                          translate=False, cpu=3, excluded_taxa=None, base_tempdir=None):
    ##TODO make the output of the translation also a dictionary? maybe?
    sequences = parse_seqs(query_fasta)
    cmds = []
    visited_queries = set()
    #print('use process_nog_hits_file')
    # if skip_queries:
    #     visited_queries.update(skip_queries)

    #tempdir = mkdtemp(prefix='emappertmp_phmmer_', dir=base_tempdir)
    ## hits_file is the output of hmmscan
    with open(parsed_hmm_output) as file:
        for line in file:
            if line.startswith('#'):
                continue

            #fields = list(map(str.strip, line.split('\t')))
            fields=line.strip().split()
            #print(fields)
            seqname = fields[0]

            if fields[1] == '-' or fields[1] == 'ERROR':
                continue

            if seqname in visited_queries:
                continue

            hitname = cleanup_og_name(fields[1])
            try:
                level = og2level[hitname]
            except KeyError:
                print(fields)
                break
            try:
                seq = sequences[seqname]
                # print(seq)
            except TypeError:
                print(seqname)
                print(type(seqname))
                break
            visited_queries.add(seqname)
            target_fasta = '/Users/zavo603/anaconda3/envs/python2/lib/python2.7/site-packages/data/OG_fasta/%s/%s.fa' % (level,hitname)
            cmds.append([seqname, seq, target_fasta, excluded_taxa])
        print('cmds',cmds)
        if cmds:
            for item in cmds:
                print(item)
                output= refine_hit(item)
                print(output)
                yield output
            # pool = multiprocessing.Pool(cpu)
            # for r in pool.imap(refine_hit, cmds):
            #     yield r
            # pool.terminate()

    #shutil.rmtree(tempdir)


def cleanup_og_name(name):
    # names in the hmm databases are sometiemes not clean eggnog OG names
    #logging.info('using cleanup_og_name')
    m = re.search('\w+\.((ENOG41|COG|KOG|arCOG)\w+)\.', name)
    if m:
        name = m.groups()[0]
    name = re.sub("^ENOG41", "", name)
    return name
###FIXME
def refine_hit(args):
    print('using refine_hits')
    seqname, seq, group_fasta, excluded_taxa = args
    #F = NamedTemporaryFile(delete=True,mode='w+t')
    F = open('refine_hit_int.txt','w')
    F.write('>%s\n%s' % (seqname, seq))
    #F.flush()
    F.close()
    best_hit = get_best_hit(F.name, group_fasta, excluded_taxa)


    return [seqname] + best_hit

## this one might need to run the shell inside the script
def get_best_hit(target_seq, target_og, excluded_taxa):
    # print(target_seq)
    # print(target_og)
    # print(excluded_taxa)
    # import subprocess
    print('using get_best_hit')
    # if not PHMMER:
    #     raise ValueError('phmmer not found in path')
    tempdir = 'get_best_hit_output'
    tempout = os.path.join(tempdir, uuid.uuid4().hex)
    #PHMMER= 'phmmer'
    cmd = "/Users/zavo603/anaconda3/bin/phmmer --incE 0.001 -E 0.001 -o rando.txt --noali --tblout %s %s %s" %\
          (tempout, target_seq, target_og)
 # file change
    # print(cmd)
    # p = subprocess.Popen("which phmmer", shell=True, stdout=subprocess.PIPE)
    # while True:
    #     line = p.stdout.readline()
    #     if line != b"":
    #         print(line)
    #     else:
    #         break
    #
    status = subprocess.call(cmd,shell=True)
    # print(status)
    best_hit_found = False
    if status == 0:
        # take the best hit
        for line in open(tempout):
            if line.startswith('#'):
                continue
            else:
                fields = line.split()
                best_hit_name = fields[0]
                best_hit_evalue = float(fields[4])
                best_hit_score = float(fields[5])
                #if not excluded_taxa or not best_hit_name.startswith("%s." % (excluded_taxa)):
                best_hit_found = True
                break
        # os.remove(tempout)
    else:
        raise ValueError('Error running PHMMER')

    if not best_hit_found:
        best_hit_evalue = '-'
        best_hit_score = '-'
        best_hit_name = '-'

    return [best_hit_name, best_hit_evalue, best_hit_score]

def main(hmmscan_output,hits_file,translated_fasta,refinement_output,og2levels):
    dump_hmm_matches(hmmscan_output, hits_file)
    refined_matches =refine_matches(translated_fasta, refinement_output, hits_file,og2levels)

if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument(
        "--hmmscan_output",
        type=str,
        help="hmmscan output file",
    )
    p.add_argument("hits_file", type=str, help="hits file output name")
    p.add_argument(
        "--translated_fasta", type=str, help="og fasta file that is translated"
    )
    p.add_argument(
        "--refinement_output",
        type=str,
        help="output file for refinement",
    )
    p.add_argument(
        "--og2levels",
        type=str,
        help="og2levels file from eggnog db",
    )

    args = p.parse_args()
    # logging.basicConfig(
    #     level=logging.INFO, datefmt="%Y-%m-%d %H:%M", format="[%(asctime)s] %(message)s"
    # )
    main(
        args.hmmscan_output,
        args.hits_file,
        args.translated_fasta,
        args.refinement_output,
        args.og2levels,
    )
