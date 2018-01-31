#-*- coding: utf - 8 - *-
#!/usr/bin/env python3

"""
Created on Monday 25 11:23:04 2017

@author: Tori,Jakob

"""

import datetime
from JsonFile import JsonFile
import json
import os
import re  # Used to retrieve directoy information
import urllib3  # Used to change HGVS with mutalyzer
import re
import csv
import pandas
import sys
import argparse
import hgvs
import hgvs.parser  # hgvs parsing and validation library
import hgvs.dataproviders.uta
import hgvs.assemblymapper
import hgvs.validator
import hgvs.exceptions
import pdb


# ===============================
# ===== helper functions ========
# ===============================


def purge(dir, pattern):
    """Deletes files in a directory which include a certain pattern"""
    for f in os.listdir(dir):
        if re.search(pattern, f):
            os.remove(os.path.join(dir, f))


def parsehgvs(hgvscode, parser, hdp):
    variant = parser.parse_hgvs_variant(hgvscode)
    vm = hgvs.assemblymapper.AssemblyMapper(
        hdp, assembly_name='GRCh37', alt_aln_method='splign')
    var_g = vm.c_to_g(variant)
    chrom = int(var_g.ac.split(".")[0][-2:])
    offset = int(var_g.posedit.pos.start.base)
    ref = var_g.posedit.edit.ref
    try:
        alt = var_g.posedit.edit.alt
    except AttributeError:
        if "dup" in str(var_g):
            alt = ref + ref
        elif "del" in str(var_g):
            alt = ""
        else:
            alt = "."
    return(chrom, offset, ref, alt)


def move(file, dir, overwrite=False):
    """
    Moves a file to a directory

    Args:
        file (string): path to file
        dir (string): path to directory
        overwrite (boolean): True if file should be replaced

    Returns:file
        exists (boolean) = False if file is already in the directory
    """
    jsonname = os.path.basename(os.path.normpath(file))
    present = os.path.isfile(dir + "/" + jsonname)
    if not present:
        os.rename(file, dir + "/" + jsonname)
    elif present and overwrite:
        os.remove(dir + "/" + jsonname)
        os.rename(file, dir + "/" + jsonname)
    return(present)


def check_hgvs(parser, validator, hdp, hgvscode):
    """Checks the HGVS code. Returns False if HGVS code is incorrect or doesnt have the correct format.

    Args:
        parser (hgvs.parser.Parser object): Used to parse the HGVS code
        validator (hgvs.validator.Validator object): Used to validate the parsed HGVS code
        hgvs (string): Input that is going to be tested

    Returns:
        correct (boolean): True if HGVS-Code is correct
    """
    try:
        variant = parser.parse_hgvs_variant(hgvscode)
        vm = hgvs.assemblymapper.AssemblyMapper(
        hdp, assembly_name='GRCh37', alt_aln_method='splign')
        var_g = vm.c_to_g(variant)
        correct = validator.validate(parser.parse_hgvs_variant(str(var_g)))
    except Exception:
        return(False)
    return(True)


def savetomultivcf(jsonobject, vcffolder, hdp):
    """
    Args:
        jsonobject (jsonobject): JsonFile object
    """

    jsondata = jsonobject.raw

    # create multivcf

    multivcf = pandas.DataFrame(columns=[
                            '#CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT', 'NM'])

    caseID = str(jsondata['case_id'])
    hgvslist = []
    rowcounter = 0
    hp = hgvs.parser.Parser()
    for mutation in jsondata['genomicData']:
        if 'HGVS-code' in mutation['Mutations']:
            hgvscode = mutation['Mutations']['HGVS-code']
            hgvslist.append(hgvscode)
        elif 'Mutation 1' in mutation['Mutations']:
            for mutationnr, mutationdict in mutation['Mutations'].items():
                if 'HGVS-code' in mutationdict:
                    hgvscode = mutationdict['HGVS-code']
                    hgvslist.append(hgvscode)
    genotype = mutation['Test Information']['Genotype']
    if genotype == 'Hemizygous':
        genotype = '1'
    elif genotype == 'Homozygous':
        genotype = '1/1'
    elif genotype == 'Heterozygous' or genotype == 'Compound Heterozygous':
        genotype = '0/1'
    else:
        genotype = './1'
    for hgvscode in set(hgvslist):
        try:
            chrom, offset, ref, alt = parsehgvs(hgvscode, hp, hdp)
            multivcf.set_value(rowcounter, '#CHROM', chrom)
            multivcf.set_value(rowcounter, 'NM', hgvscode)
            multivcf.set_value(rowcounter, 'POS', offset)
            multivcf.set_value(rowcounter, 'ID', '.')
            multivcf.set_value(rowcounter, 'REF', ref)
            multivcf.set_value(rowcounter, 'ALT', alt)
            multivcf.set_value(rowcounter, 'QUAL', '.')
            multivcf.set_value(rowcounter, 'FILTER', '.')
            multivcf.set_value(rowcounter, 'INFO', 'HGVS="' + hgvscode + '"')
            multivcf.set_value(rowcounter, 'FORMAT', 'GT')
            multivcf.set_value(rowcounter, caseID, genotype)
            rowcounter = rowcounter + 1
        except ValueError:
            print ('HGVS-Code not parsed:', caseID, hgvscode)
            continue

    # data_vcf sortieren
    multivcf = multivcf.sort_values(by=['#CHROM', "POS"])
    multivcf = multivcf.reset_index(drop=True)
    multivcf = multivcf.fillna(value='0/0')

    with open(vcffolder + "/" + caseID + ".vcf", 'a') as outfile:
        outfile.write(
            '##fileformat=VCFv4.1\n##INFO=<ID=HGVS,Number=1,Type=String,Description="HGVS-Code">\n##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n')

    multivcf.to_csv(vcffolder + "/" + caseID + ".vcf", mode="a",sep='\t', index=False,
                    header=True, quoting=csv.QUOTE_NONE)


# ===============================
# ===== Inspection function =====
# ===============================

def checkjson(jsonobject, errordict, hdp):
    """Checks an jsonobject for missing values and the HGVS-Code.

    Args:
        jsonobject(JsonFile): A JsonFile object


    Returns:
        correct (boolean): True if json-file is correct
        errordesc (list): error descriptions
        jsonobject (JsonFile): The input JsonFile object with corrected HGVS-Codes (if errorcorrection was applied)
    """

    errordesc = []
    correct = True

    with open(errordict) as json_data:
        transcript_errors = json.load(json_data)

    d = jsonobject.raw

    vcfexists = d['vcf'] != 'noVCF'

    if isinstance(d['geneList'], list):
        gscore_list = []
        for gene in d['geneList']:
            if "gestalt_score" in gene:
                gscore = gene['gestalt_score']
                gscore_list.append(float(gscore))
        if gscore_list:
            max_gscore = max(gscore_list)
    else:
        correct = False
        errordesc.append("no geneList")

    # Check for annotated features
    if not d['features']:
        correct = False
        errordesc.append("no features")

    # No image uploaded
    if max_gscore == 0:
        correct = False
        errordesc.append("no Picture")

    # Check for the molecular diagnosis
    if not d["ranks"]:
        correct = False
        errordesc.append("no diagnosis")

    # HGVS-code check

    # Initialize HGVS parser and validator
    parser = hgvs.parser.Parser()
    validator = hgvs.validator.Validator(hdp=hdp)

    if not d['genomicData']:
        correct = False
        errordesc.append("no genomicData")
    for idx1, mutation in enumerate(d['genomicData']):
        if not mutation['Mutations']:
            errordesc.append("one of the mutations not listed")
            continue
        if 'HGVS-code' in mutation['Mutations']:
            hgvscode = mutation['Mutations']['HGVS-code']
            hgvscheck = check_hgvs(parser, validator, hdp, hgvscode)
            if not(hgvscheck):
                if hgvscode in transcript_errors:
                    if(check_hgvs(parser,validator,hdp,transcript_errors[hgvscode])):
                        d["genomicData"][idx1]["Mutations"]["HGVS-code"] = transcript_errors[hgvscode]
                        continue
                    else:
                        correct = False
                        errordesc.append("Wrong HGVS-Code: " + hgvscode)
                else:
                    correct = False
                    errordesc.append("Wrong HGVS-Code: " + hgvscode)

        # Check for compound heterozygous
        elif 'Mutation 1' in mutation['Mutations']:
            for idx2, mutationdict in enumerate(mutation['Mutations']):
                if 'HGVS-code' in mutationdict:
                    hgvscode = mutationdict['HGVS-code']
                    hgvscheck, errorstring = check_hgvs(
                        parser, validator, hdp, hgvscode)
                    if not(hgvscheck):
                        if hgvscode in transcript_errors:
                            if(check_hgvs(parser,validator,hdp,transcript_errors[hgvscode])):
                                d["genomicData"][idx1]["Mutations"][idx2]["HGVS-code"] = transcript_errors[hgvscode]
                                continue
                            else:
                                correct = False
                                errordesc.append("Wrong HGVS-Code: " + hgvscode)
                        else:
                            correct = False
                            errordesc.append(
                                "Wrong HGVS-Code: " + hgvscode)

        # Check for monogenic
        elif mutation['Test Information']['Mutation Type'] == 'Monogenic':
            if 'protein level' in mutation['Mutations']:
                correct = False
                errordesc.append("Mutation on protein level")
            else:
                if not vcfexists:
                    correct = False
                    errordesc.append("Monogenic and no vcf file listed")
                else:
                    correct = False
                    errordesc.append("VCF present but no HGVS-Code")
        else:
            correct = False
            errordesc.append("No HGVS-Code listed")
        if mutation['Test Information']['Mutation Type'] != 'Monogenic':
            correct = False
            errordesc.append("Not Monogenic")
    return(correct, errordesc, jsonobject)


# ===============================
# == errorcorrection functions===
# ===============================
    #
    # def errorcorrbymut(hgvscode,jsonfile): #TODO for what files do we need this?
    #     """
    #     Tries to correct hgvs code without transcript with the mutalyzer
    #
    #     Args:
    #         hgvscode (string): hgvscode with missing transcript
    #
    #     Returns:
    #         correct(boolean): True if json file was corrected
    #         hgvscode(string): containing either the old hgvscode (if not corrected) or the corrected version
    #     """
    #
    #     with open(jsonfile) as json_file:
    #         jsondata = json.load(json_file)
    #         url = "https://mutalyzer.nl/position-converter?assembly_name_or_alias=GRCh37&description=" + hgvscode
    #         try:
    #             page = urllib2.urlopen(url)
    #             data = page.read()
    #         except Exception:
    #             return(False)
    #         if 'Found transcripts' in data:
    #             data = data.split(
    #                 'variant region</h4>')[1].split('<br></pre>')[0]
    #             transcript = key.split('.')[0]
    #             transcriptversion = key.split(
    #                 ':')[0].split('.')[1]
    #             positions = [m.start()
    #                          for m in re.finditer(transcript, data)]
    #             alt_transcripts = []
    #             for position in positions:
    #                 trans = data[position:].split(':')[0]
    #                 trans_vers = trans.split(transcript)[1].split(':')[
    #                     0].split('.')[1]
    #                 if trans_vers > transcriptversion:
    #                     transcriptversion = trans_vers
    #             transcriptversion = transcript + '.' + \
    #                 transcriptversion + ':' + key.split(':')[1]
    #             for mutation in d['genomicData']:
    #                 if 'HGVS-code' in mutation['Mutations']:
    #                     if mutation['Mutations']['HGVS-code'] == hgvscode:
    #                         mutation['Mutations']['HGVS-code'] = transcriptversion
    #                 elif 'Mutation 1' in mutation['Mutations'].keys():
    #                     for multimut, description in mutation['Mutations'].items():
    #                         if mutation['Mutations'][multimut]['HGVS-code'] == hgvscode:
    #                             mutation['Mutations']['HGVS-code'] = transcriptversion
    #             json.dump(jsondata,json_file)
    #             checkjson(jsonfile)
    #         elif 'We found these versions' in data:
    #             # print jsonfile, data
    #             newtranscript = data.split('We found these versions: ')[
    #                 1].split('<p></p>')[0].split('</p>')[0]
    #             newtranscript = newtranscript + \
    #                 ':' + key.split(':')[1]
    #             with open(jsoncurratedfolder + '/' + jsonfile) as json_data:
    #                 d = json.load(json_data)
    #                 for mutation in d['genomicData']:
    #                     if 'HGVS-code' in mutation['Mutations'].keys():
    #                         if mutation['Mutations']['HGVS-code'] == key:
    #                             mutation['Mutations']['HGVS-code'] = newtranscript
    #                             # append_errordict()
    #                         elif 'Mutation 1' in mutation['Mutations'].keys():
    #                             for multimut, description in mutation['Mutations'].items():
    #                                 if mutation['Mutations'][multimut]['HGVS-code'] == key:
    #                                     mutation['Mutations']['HGVS-code'] = newtranscript
    #                                     # append_errordict()
    #             transcript_errors[key] = newtranscript
    #             muterrors = muterrors + 1
    #             print jsonfile,  'vorher: ', key, 'nachher: ', newtranscript
    #             new_hgvs.append(newtranscript)
    #             corrected_transcripts[key] = newtranscript
    #         elif 'could not be found in our database (or is not a transcript).' in data:
    #             print 'no transcript found: ', submitter, jsonfile, key
    #         else:
    #             print 'please check: ', jsonfile, key
    #
    #         with open(errordict, 'w') as dicttojson:
    #             json.dump(transcript_errors, dicttojson)
    #
    #         step = 'mutalyzer'
    #         overview = checkjsons(step, jsoncurratedfolder)
    #
    # def errorcorrmanual(jsonfile, errordict, hdp):
    #     """
    #     Tries to correct errors manually
    #
    #     Args:
    #         jsonfile (string): path to jsonfile
    #         errordict (string): path to errordictionary
    #
    #     Returns:
    #         correct (boolean): True if the json file has been corrected
    #         error (string): error description
    #     """
    #     with open(jsonfile, 'r') as jsonfile:
    #         jsondata=json.load(jsonfile)
    #         for mutation in jsondata['genomicData']:
    #             if 'HGVS-code' in mutation['Mutations']:
    #                 print(mutation['Mutations']['HGVS-code'])
    #                 right_hgvs=input('right HGVS (u for unknown):')
    #                 if right_hgvs == 'u':
    #                     continue
    #                 else:
    #                     mutation['Mutations']['HGVS-code']=right_hgvs
    #                     transcript_errors[hgvs]=right_hgvs
    #                     with open(errordict, 'w') as dicttojson:
    #                         json.dump(transcript_errors, dicttojson)
    #             elif 'Mutation 1' in mutation['Mutations']:
    #                 for multimut, description in mutation['Mutations'].items():
    #                     print(mutation['Mutations'][multimut]['HGVS-code'])
    #                     right_hgvs=input('right HGVS (u for unknown):')
    #                     if right_hgvs == 'u':
    #                         continue
    #                     else:
    #                         mutation['Mutations'][multimut]['HGVS-code']=right_hgvs
    #                         transcript_errors[hgvs]=right_hgvs
    #                     with open(errordict, 'w') as dicttojson:
    #                         json.dump(transcript_errors, dicttojson)
    #             with open(jsonfile, 'w') as dicttojson:
    #                 json.dump(jsondata, dicttojson)
    #             correct, error=checkjson(jsonfile, hdp)
    # return(correct, error)


# ===============================
# ===== main script =============
# ===============================

def main():
    """
    main function with the following steps:
    1. Read all inputs
    2. Initialize JsonFile object of the input jsonfile
    3. Check for errors
        3a. Create vcf for json file and move file to curratedfolder if no errors occured
        3b. Move jsonfile to the debugfolder and add error descriptions to json dict containing the error description
    """

    parser = argparse.ArgumentParser(description='jsonfile-qualitycheck')
    parser.add_argument('--jsonfile', metavar='pathtojsonfile', type=str, nargs="?",
                        help='Path to the json-file')
    parser.add_argument('--jsoncurratedfolder', metavar='pathtojsoncuratedfolder', type=str, nargs="?",
                        help='Path to the folder containing the currated json-files')
    parser.add_argument('--errordict', metavar='pathtoerrordict', type=str, nargs="?",
                        help='Path to the HGVS-errordictionary in json-format')
    parser.add_argument('--logfile', metavar='pathtologfile', type=str, nargs="?",
                        help='Path to the logfile')
    parser.add_argument('--debugfolder', metavar='pathtodebugfolder', type=str, nargs="?",
                        help='Path to the Folder containing Errorreports')
    parser.add_argument('--vcffolder', metavar='pathtovcffolder', type=str, nargs="?",
                        help='Path to the v')

    args = parser.parse_args()

    jsonfile = args.jsonfile
    jsoncurratedfolder = args.jsoncurratedfolder
    errordict = args.errordict
    logfile = args.logfile
    debugfolder = args.debugfolder
    vcffolder = args.vcffolder

    # Initialize JsoFile object
    jsonobject = JsonFile(jsonfile)

    # Initialize server for HGVS-code validation
    hdp = hgvs.dataproviders.uta.connect()

    # check if json-file is correct
    correct, errordesc, jsonobject = checkjson(jsonobject, errordict, hdp)

    now = datetime.datetime.now()
    date = now.strftime('%Y-%m-%d')
    time = now.strftime('%H:%M:%S')

    if correct:  # if the file has passed save it to to a vcf file and move it to the currated folder
        savetomultivcf(jsonobject, vcffolder, hdp)
        with open(jsonfile,"w") as newjson:
            json.dump(jsonobject.raw,newjson)
        present = move(jsonfile, jsoncurratedfolder)
        if present:
            overwrite = input(
                "This jsonfile has already been currated \n overwrite?(y or n):")
            if overwrite == "y":
                move(jsonfile, jsoncurratedfolder, overwrite=True)
    else:  # if errors occured move file to the debug folder and add errors to the errordictionary in the debug folder
        with open(debugfolder + "/errorsummary.json", "r+") as errorlist:
            errordata=json.load(errorlist)
            errordata.update({jsonobject.raw["case_id"]: errordesc})
            errorlist.seek(0)
            json.dump(errordata,errorlist)
            errorlist.truncate()
        present = move(jsonfile, debugfolder)
        if present:
            overwrite = input(
                "This jsonfile is already in the debug fodler \n overwrite?(y or n):")
            if overwrite == "y":
                move(jsonfile, jsoncurratedfolder, overwrite=True)

    # save to logfile
    progress = pandas.read_excel(open(logfile, 'rb'), sheetname=0)
    newest = [date, time, jsonobject.raw["case_id"], str(correct)]
    progress.loc[len(progress)] = newest
    progress.to_excel(logfile, index=False)


if __name__ == "__main__":
    main()
