import svist4get as sv4g
import pandas as pd
import re
import os
from IPython.display import Image, display
#set up the main params
def init_pa(exp_dict,exp='RNAseqInVivo_dh'):
    pa = sv4g.manager.Parameters()
    path_to_config = exp_dict[exp]['path_to_config']
    pa.initialize(path_to_config)
    pa.config['gtf_file'] = exp_dict[exp]['gtf_file']
    pa.config['fasta_file'] = exp_dict[exp]['fasta_file']
    pa.config['bedgraph'] =  exp_dict[exp]['paths_to_bedgraphs']
    pa.config['bedgraph_label'] = exp_dict[exp]['bedgraph_label']
    pa.config['bedgraph_label_position'] = 'left'
    pa.config['transcript_label_style'] = 'gene_id'
    #pa.config['bedgraph_upper_limit'] = 'max'
    return pa

#parse the gtf file so that can be serched
def make_info_dict(info,sep='; '):
    def replace_chr(line):
        line = re.sub('[;\"]', '', line)
        return line
    
    res = {}
    for item in info.split(sep):
        #print(item.split(' '))
        if len(item.split(' '))>1:
            res[ replace_chr(item.split(' ')[0]) ]=replace_chr(item.split(' ')[1].strip())
    return res

def parse_gtf(pa):
    gtf = pd.read_csv(pa.config['gtf_file'], sep='\t', header=None,comment='#')
    gtf.columns = ['chro','source','ftype','start','end','score','strand','frame','info']
    #print(gtf.head())
    #gtf['gene'] = [n.split(';')[0].split(' ')[-1].split(':')[0].strip('\"') for n in gtf['info']]
    #temp = []
    #for info in gtf['info']:
        #print(info)
        #make_info_dict(info)['gene_id']
    gtf['gene'] = [make_info_dict(info)['gene_id'] for info in gtf['info']]
    #print(gtf.tail())
    return gtf

#extract the gene of interest +- n genes  
#at 5' and 3' of the gene of interest
def get_region(in_gene, pa, extend=2):
    gtf = parse_gtf(pa)
    temp = gtf.drop_duplicates(subset=['chro','gene'])
    temp = temp[temp.iloc[:,2]=='transcript']
    #print(temp[temp['gene'].str.contains(in_gene)])
    #print( temp[temp['gene']==in_gene])
    selection = temp[temp['gene']==in_gene].index.values[0]
    
    if extend == 0:
        chrom = temp.loc[selection]['chro']
        start=temp.loc[selection]['start']
        end =temp.loc[selection]['end']
        strand=temp.loc[selection]['strand']
        #print(chrom, start, end, strand)
        #print('im here 3')
        return chrom, start, end, strand
   
    #print(selection)
    #print('____________')
    temp = temp[temp['chro']==temp.loc[selection]['chro']]
    temp = temp.sort_values('start')
    
    #print(temp.head())
    #print('____________')
    #print(temp.tail())
    #print('____________')
    
    strand = temp.loc[selection]['strand']
    chrom = temp.loc[selection]['chro']
    temp['old_index']=temp.index.values
 
    
    temp = temp.reset_index(drop=True)

    #print(temp.head())
    #print('____________')       
    
    new_pos = temp[temp.old_index==selection].index.values[0]
    
    
    from_index = new_pos-extend-1
    to_index = new_pos+extend+1
    #print(temp.head())
    #print('____________')
    #print(temp.loc[selection])
    
    #we ask for an array of gene
    #we take n genes at the 5' and n genes at the 3' of the selected gene where n=extend 
    #we also add one further gene at the start and end of the array
    #print(from_index, to_index)
    
    #check if we are at the start of the chr
    if from_index<0:
        #print('-------')
        #print(from_index)
        from_index=0
        temp = temp.loc[from_index:to_index]
        #start of plot is the start of the gene of interest
        #not ideal
        start = temp['start'].values[0]
        end = temp['start'].values[-1]
        #print(start,end)
        #print('im here 5')
        return chrom, start, end, strand
    
    #or end of the chromosome    
    elif to_index > temp.index.values[-1]:
        #print('-------')
        
        #print(to_index, temp.index.values[-1])
        to_index = temp.index.values[-1]
        temp = temp.loc[from_index:to_index]
        #print(temp)
        
        start = temp['end'].values[0]
        #end of plot is the end of the gene of interest
        #not ideal        
        end = temp['end'].values[-1]
        #print(start,end)
        #print('im here 4')
        return chrom, start, end, strand

    else:
        #otherwise    
        temp = temp.loc[from_index:to_index]
        #the start of the plot is the end of the first gene in the array
        #start = temp['end'].values[0]
        start = temp['start'].values[1]
        #the end of the plot is the start of the last gene in the array
        #end = temp['start'].values[-1]
        end = temp['end'].values[-2]
        #print(temp)
        #print(1, 'extend', extend)
        if start > end:
            #print('im here 1')
            #print(chrom, end, start, strand)
            return chrom, end, start, strand
        #print('im here 2')
        #print(chrom, start, end, strand)    
        return chrom, start, end, strand

#set up the gene-dependent parameters
def add_gene(gene, desc,  pa, extend=2, out='', always_forward=False):
    chrom, start, end, strand = get_region(gene, pa, extend=extend) 
    #print(chrom, start, end, strand)
    pa.config['window'] = [chrom, start, end]
    pa.config['image_title'] = gene+' '+desc
    pa.config['output_filename'] = os.path.join(out, gene)
    #determine the orientation of the plot
    if always_forward:
        if strand == '+':
            pa.config['revcomp_transform'] = 0
        else:
            pa.config['revcomp_transform'] = 1
    
    return pa
    #print(2, 'extend', extend)

def add_reagion(chrom, start, end, title, pa, out):
    pa.config['window'] = [chrom, start, end]
    pa.config['image_title'] = title
    pa.config['output_filename'] = os.path.join(out, gene)
    return pa

# make the figure
def make_image(pa, do_display=False, add_aa=False, add_transcripts=True):
    #print("pa.config['gtf_file']", pa.config['gtf_file'])
    gtf = sv4g.data_processing.Gtf_helper(pa.config['gtf_file'])
    transcripts = gtf.extract_transcripts_from_widnow(*pa.config['window'])
    data_from_gtf = (gtf.extract_data_about_transcripts(transcripts))
    #print(data_from_gtf)
    pa.add_gtf_data(data_from_gtf)
    tracks = []
    tracks += sv4g.manager.Title_tracks_maker(pa).create_tracks()
    tracks += sv4g.manager.Axis_tics_tracks_maker(pa).create_tracks()
    tracks += sv4g.manager.Vgrid_tracks_maker(pa).create_tracks()
    if add_aa:
        tracks += sv4g.manager.Aa_seq_tracks_maker(pa).create_tracks()
    if add_transcripts:
        tracks += sv4g.manager.Transcript_struct_tracks_maker(pa).create_tracks()
    tracks += sv4g.manager.Bedgraph_tracks_maker(pa).create_tracks()
    sv4g.manager.Image(tracks, pa).draw()
    # converting the resulting pdf to a png file
    sv4g.methods.pdf_page_to_png(pa)
    if do_display:
        display(Image(filename=os.path.join(pa.config['output_filename']+'.png')))

def find_max_coverage_in_gene(gene, pa, bg_df_list):
    #get gene coordinates
    #print(1)
    region = get_region(gene, pa, extend=0)
    #select from each bg dataframe the reagion of interest
    selection = [ df[(df['start']>region[1]) & (df['end']<region[2]) & (df['chr']==region[0])] for df in bg_df_list]
    #concat for all bg trak and take max coverage
    selection = pd.concat(selection)
    #print(2)
    return selection['cov'].max()

def find_max_coverage_in_region(gene, pa, bg_df_list,extend=0):
    #get gene coordinates
    #print(1)
    region = get_region(gene, pa, extend=extend)
    #select from each bg dataframe the reagion of interest
    selection = [ df[(df['start']>region[1]) & (df['end']<region[2]) & (df['chr']==region[0])] for df in bg_df_list]
    #concat for all bg trak and take max coverage
    selection = pd.concat(selection)
    #print(2)
    return selection['cov'].max()
