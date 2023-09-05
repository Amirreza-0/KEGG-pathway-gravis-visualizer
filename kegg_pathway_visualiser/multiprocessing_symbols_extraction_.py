import json
import os
import numpy as np
import requests
import multiprocessing
from joblib import Parallel, delayed
import re
from Bio.KEGG.KGML.KGML_parser import read

def extract_symbols(pathway_file_name):
    '''
    This function is used to extract the symbols from the KEGG website. The symbols are used to replace the gene IDs
    in the pathway. The symbols are stored in a json file to avoid extracting them again.
    param pathway_file_name: the name of the pathway file
    :return: the path to the symbols json file
    '''

    def load_previous_symbols():
        if not os.path.exists(f'media/symbols.json'):
            with open(f'media/symbols.json', 'w') as fp:
                json.dump({}, fp)
        with open(f'media/symbols.json', 'r') as fp:
            gene_IDs = json.load(fp)
        return gene_IDs

    prev_symbols = load_previous_symbols()

    def load_pathway(pathway_name):
        pathway = read(open(pathway_name, 'r'))
        return pathway

    pathway = load_pathway(f'media/{pathway_file_name}')


    def extract_symbol(link, indx=0):
        # extracting the sybmol from the link attribute of the gene
        doc = requests.get(link).text
        # everything between <td class="td11 defd"><div class="cel"> and <br> </div></td>
        pattern = re.compile(r'<td class="td11 defd"><div class="cel">(.*)<')
        # find all the matches
        symbols = pattern.findall(doc)[0]
        # choosing the symbol by index (since there could be multiple symbols in one link)
        symbol = symbols.split(", ")[indx]
        return symbol


    def process_gene_ids(gene_ids, link):
        final_dict = {}
        for indx, g in enumerate(gene_ids):
            human_readable_name = extract_symbol(link, indx)
            final_dict[g] = human_readable_name
            #print(f"original_name: {g}, human_readable_name: {human_readable_name}")
        return final_dict


    def dig(gene, final_dict):
        print(f"working on link: {gene.link}")

        try:
            gene_id = gene.name.split(" ")
            for g in gene_id:
                if g in prev_symbols.keys():
                    final_dict[g] = prev_symbols[g]
                    return final_dict

        except:
            gene_id = gene.name
            if gene_id in prev_symbols.keys():
                final_dict[gene_id] = prev_symbols[gene_id]
                return final_dict
            pass


        try:
            if len(gene_id) > 10:
                new_links = []
                ids = gene_id
                while len(ids) > 10:
                    temp = ids[:10]
                    ids = ids[10:]

                    new_links.append(f"https://www.kegg.jp/dbget-bin/www_bget?{'+'.join(temp)}")
                new_links.append(f"https://www.kegg.jp/dbget-bin/www_bget?{'+'.join(ids)}")

            if len(gene_id) > 1:
                #print(f"gene family found with length {len(gene_id)}")
                if len(gene_id) < 10:
                    final_dict.update(process_gene_ids(gene_id, gene.link))
                if len(gene_id) > 10:
                    for link in new_links:
                        while len(gene_id) > 10:
                            temp = gene_id[:10]
                            gene_id = gene_id[10:]
                            final_dict.update(process_gene_ids(temp, link))
                        try:
                            final_dict.update(process_gene_ids(gene_id, link))
                        except:
                            pass

            human_readable_name = extract_symbol(gene.link)
            final_dict[gene_id[0]] = human_readable_name
            #print(f"original_name: {gene_id}, human_readable_name: {human_readable_name}")

        except:
            pass

        return final_dict


    final_dict = {}

    items = []
    for entry_id in pathway.entries:
        if pathway.entries[entry_id].name not in prev_symbols.keys():
            if pathway.entries[entry_id].name != "undefined":
                items.append(pathway.entries[entry_id])

    for relation in pathway.relations:
        if relation.entry1.name not in prev_symbols.keys():
            items.append(relation.entry1)
        if relation.entry2.name not in prev_symbols.keys():
            items.append(relation.entry2)

    print(len(prev_symbols), "symbols already extracted")
    print("extracting symbols for: ", len(items), "new items")

    # applying multiprocessing to speed up the process
    number_of_cores = multiprocessing.cpu_count()
    with Parallel(n_jobs=number_of_cores) as parallel:
        output = parallel(delayed(dig)(itm, final_dict) for itm in items)
        outputs = np.array(output)
        for out in outputs:
            final_dict.update(out)

    prev_symbols.update(final_dict)

    with open(f'media/symbols.json', 'w') as fp:
        json.dump(prev_symbols, fp)

    symbols_path = f'media/symbols.json'

    return symbols_path

