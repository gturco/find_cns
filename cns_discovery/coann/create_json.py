import json

def create_json(query,subject,blast,outdir):
  "formats my file information into a json file needed for coannot"
  json_file = open("{2}/{0}_{1}.json".format(query,subject,outdir), "w")
  json_dict = {"genome_a":{"name":"{0}".format(query),"flat":"{1}/{0}.bed".format(query,outdir),
  "fasta":"{1}/{0}.fasta".format(query,outdir)},"genome_b":{"name":"{0}".format(subject),"flat":"{0}/{1}.bed".format(outdir,subject),
  "fasta":"{0}/{1}.fasta".format(outdir,subject)},"blast": {"W": 20,"a": 8,"e": 0.001},
  "default":{"blast_path":"{0}/blastall".format(blast),"out_dir":"{0}".format(outdir),"min_len":100, "reciprocal":True,"blast_log":True}}
  json_formatted = json.dumps(json_dict,indent=4)
  json_file.write(json_formatted)
  json_file.close()

if __name__ == "__main__":
    import optparse
    parser = optparse.OptionParser()
    parser.add_option("--query", dest="query", help="name of query org or ORGA")
    parser.add_option("--subject", dest="subject", help="name of subject org or ORGB")
    parser.add_option("--blast_path",dest="blast", help="location of most up todate blast dir")    
    parser.add_option("--out_dir", dest="outdir", help="dir were data is kept")
    (options, _) = parser.parse_args()
    create_json(options.query, options.subject,options.blast,options.outdir)

