import pandas as pd
import numpy as np

import gspread
import glob as glob
import re 

import statsmodels.formula.api as smf
import statsmodels.api as sm
import statsmodels.sandbox.stats.multicomp as mc

import matplotlib.pyplot as plt
from matplotlib import cm, colors, patches, ticker
from matplotlib.figure import Figure
from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas
import seaborn as sns

import json
from flask import jsonify

try:
    from StringIO import StringIO
except ImportError:
    from io import StringIO
 
import markdown

from flask import Flask, make_response, render_template, Markup, request, redirect, Response
#from flask_bootstrap import Bootstrap

from werkzeug.contrib.fixers import ProxyFix

import requests

#import urllib2, urllib
try:
    from urllib.request import urlopen, HTTPCookieProcessor
except ImportError:
    # Fall back to Python 2's urllib2
    from urllib2 import urlopen, HTTPCookieProcessor

#try:
#    # For Python 3.0 and later
#    from http.cookiejar import CookieJar
#except ImportError:
#    from cookielib import CookieJar
# import poster

# from flask_sitemap import Sitemap

app = Flask(__name__, static_url_path='/static')

# Bootstrap(app)
# ext = Sitemap()
# from flask_shorturl import ShortUrl
# surl = ShortUrl(app)
#from flask.ext.autodoc import Autodoc
#auto = Autodoc(app)


@app.route('/test_d3/<gene_of_interest>')
def return_test_d3(gene_of_interest):
    this_title = gene_of_interest
    return render_template('test_d3_gene.html', **locals())


@app.route('/parse_features_for_d3/<gene_of_interest>.csv')
def parse_features_for_d3(gene_of_interest):
    column_regex=','+gene_of_interest+'@'
    # print column_regex
    these_results = results
    these_results.index = these_results['name']
    these_results = these_results.transpose()
    these_results = these_results.filter(regex=column_regex )

    some_results = these_results.transpose().filter(regex='^Anxious Temperament \(mean\)')

    def generate():
        yield 't,p,chromasome,gene_id,gene_symbol,feature_type,annotation_type,start,stop,start2,stop2\n'
        for feature in d.filter(regex=column_regex).columns:
            # amt = these_results.ix[these_results.index==feature]['average_exprs'].values[0]
            t_value = np.round(some_results.ix[some_results.index==feature].values[0,0], decimals=2)
            p_value = np.round(some_results.ix[some_results.index==feature].values[0,1], decimals=4)

            this_data = [str(t_value),str(p_value)]

            p = re.split(r'[,@\:\=\-\_]', feature)   
    
            gene_name = p[2]
            feature_type = p[3]
            coords = [int(coord) for coord in p[5:len(p)] ]
            start = 1.0*coords[0]
            end = 1.0*coords[-1]
            middle = start+(end-start)/2.0
    
            yield ','.join(this_data) + ',' + ','.join(p) + '\n'

    return Response(generate(), mimetype='text/csv')


@app.route('/parse_features_for_d3/<gene_of_interest>.json')
def parse_features_for_d3_json(gene_of_interest):
    column_regex=','+gene_of_interest+'@'
    # print column_regex
    these_results = results
    these_results.index = these_results['name']
    these_results = these_results.transpose()
    these_results = these_results.filter(regex=column_regex )

    some_results = these_results.transpose().filter(regex='^Anxious Temperament \(mean\)')

    print( [re.split(r'[,@\:\=\-\_]', feature)[0:5][:] for feature in some_results.index] )


    # [some_results['chromasome'],some_results['gene_id'],some_results['gene_symbol'],some_results['feature_type'], some_results['annotation_type']] = [re.split(r'[,@\:\=\-\_]', feature)[0:5] for feature in some_results.index] 
    # (chromasome,gene_id,gene_symbol,feature_type,annotation_type) = [re.split(r'[,@\:\=\-\_]', feature)[0:5] for feature in some_results.index] 
    arr = [re.split(r'[,@\:\=\-\_]', feature)[0:5] for feature in some_results.index] 
    print( arr[1] )

    # for feature in d.filter(regex=column_regex).columns:
    #     # amt = these_results.ix[these_results.index==feature]['average_exprs'].values[0]
    #     some_results['t'] = np.round(some_results.ix[some_results.index==feature].values[0,0], decimals=2)
    #     some_results['p'] = np.round(some_results.ix[some_results.index==feature].values[0,1], decimals=4)


    #     p = re.split(r'[,@\:\=\-\_]', feature)   
    
    #     some_results['chromasome'] = p[0]
    #     some_results['gene_id'] = p[1]
    #     some_results['gene_symbol'] = p[2]
    #     some_results['feature_type'] = p[3]
    #     some_results['annotation_type'] = p[4]
    #     coords = [int(coord) for coord in p[5:len(p)] ]
    #     some_results['coords'] = ','.join(p[5:len(p)])
    #     some_results['start'] = 1.0*coords[0]
    #     some_results['stop'] = 1.0*coords[-1]
    #     some_results['middle'] = coords[0]+(coords[-1]-coords[0])/2.0


    print(  some_results )
    return Response(some_results.to_json(orient='index'), mimetype='application/json')




@app.route('/')
def return_welcome():
    this_title = 'Fox et al., RNA-seq and Anxious Temperament (AT)'
    option_list = [c[:-2] for c in results.filter(regex='_p$').columns]
    return render_template('index.html', **locals())

@app.route('/index.html')
def also_welcome():
    this_title = 'Kalin-Knowles RNAseq'
    option_list = [c[:-2] for c in results.filter(regex='_p$').columns]
    return render_template('index.html', **locals())

@app.route('/rnaseq')
def also_also_welcome():
    this_title = 'Kalin-Knowles RNAseq'
    option_list = [c[:-2] for c in results.filter(regex='_p$').columns]
    return render_template('index.html', **locals())

@app.route('/error', methods=['GET'])
def error_welcome(error_text='Ooops! Something has gone wrong, please try again.'):
    search_text = request.args.get('error_text', default = '', type = str)
    if search_text != '':
        error_text = error_text+'\nCannot find gene \"'+search_text+'\"'
    this_title = 'Kalin-Knowles RNAseq'
    option_list = [c[:-2] for c in results.filter(regex='_p$').columns]
    return render_template('index.html', **locals())

@app.route('/about')
def about_this_site():
    this_title = 'Kalin-Knowles RNAseq'
 
    return render_template('about_this_site.html', **locals())


@app.route('/search')
##@auto.doc()
def search_for_top():
    this_title = 'Search for Top Associations'
    option_list = [c[:-2] for c in results.filter(regex='_p$').columns]
    return render_template('search.html', **locals())

@app.route('/search_top', methods=['GET'])
def search_top(genes_or_features='features', selected_column='T1T3', expression_threshold=0, sort_col=1, n=10, threshold='.1'):
    url = '/top_'
    if 'threshold' in request.args:
        threshold = request.args.get('threshold')
    if 'genes_or_features' in request.args:
        genes_or_features = request.args.get('genes_or_features')
    if 'top_n' in request.args:
        n = request.args.get('top_n')
    if 'min_exprs' in request.args:
        expression_threshold = request.args.get('min_exprs')
    if 'selected_column' in request.args:
        selected_column = request.args.get('selected_column')
    if 'sort_col' in request.args:
        sort_col = request.args.get('sort_col')
    else:
        print(genes_or_features)
        if genes_or_features=='genes':
            print('this worked')
            sort_col=1
        else:
            sort_col=2

    url = url+genes_or_features+'/'+selected_column+'?n='+str(n)+'&min_exprs='+str(expression_threshold)+'&sort_col='+str(sort_col)+'&threshold='+str(threshold)

    return redirect(url)



@app.route('/documentation')
def documentation():
#    return auto.html()
    return auto.html(template='docs.html')


@app.route('/david_list/<gene_id_list>', methods=['GET'])
#@auto.doc()
def david_list(gene_id_list):
    if 'analysis_name' in request.args:
        analysis_name = request.args.get('analysis_name')

    shared_url_head = "http://david.abcc.ncifcrf.gov/api.jsp?type=ENTREZ_GENE_ID&ids="
    unique_url = gene_id_list
    shared_url_tail = "&tool=summary"

    davidLink = shared_url_head+unique_url+shared_url_tail

    return redirect(davidLink)



@app.route('/enrichr_list/<gene_list>', methods=['GET'])
#@auto.doc()
def enrichr_list(gene_list):
    if 'analysis_name' in request.args:
        analysis_name = request.args.get('analysis_name')
    else:
        analysis_name = 'abba'

    ENRICHR_URL = 'http://amp.pharm.mssm.edu/Enrichr/addList'
    genes_str = '\n'.join( gene_list.split(';') )
    payload = {
        'list': (None, genes_str),
        'description': (None, analysis_name)
    }

    response = requests.post(ENRICHR_URL, files=payload)
    if not response.ok:
        raise Exception('Error analyzing gene list')
    
    data = json.loads(response.text)

    shareUrlHead = "http://amp.pharm.mssm.edu/Enrichr/enrich?dataset="
    enrichrLink = shareUrlHead + data['shortId']
    return redirect(enrichrLink)



    print(data)


# def enrichr_list_depreciated(gene_list):
#     print('what?')
#     if 'analysis_name' in request.args:
#         analysis_name = request.args.get('analysis_name')
#     else:
#         analysis_name = 'abba'
# 
#     genesStr = '\n'.join( gene_list.split(';') )
#     #post a gene list to enrichr server and get the link.
#     cj = CookieJar()
#     opener = poster.streaminghttp.register_openers()
#     opener.add_handler(HTTPCookieProcessor(CookieJar()))
# 
#     params = {'list':genesStr,'description':analysis_name}
#     datagen, headers = poster.encode.multipart_encode(params)
#     url = "http://amp.pharm.mssm.edu/Enrichr/enrich"
#     enrichr_request = Request(url, datagen,headers)
#     urlopen(enrichr_request)
# 
#     x = urlopen("http://amp.pharm.mssm.edu/Enrichr/share")
#     responseStr = x.read()
#     splitPhrases = responseStr.split('"')
#     linkID = splitPhrases[3]
#     shareUrlHead = "http://amp.pharm.mssm.edu/Enrichr/enrich?dataset="
#     enrichrLink = shareUrlHead + linkID
#     cj.clear_session_cookies()
#     return redirect(enrichrLink)

def get_top_list_args( this_args ):
    expression_threshold=0
    sort_col = 2
    thr = None
    n = None
    if 'n' in this_args:
        n = int(this_args.get('n'))
    if 'min_exprs' in this_args:
        expression_threshold = int(this_args.get('min_exprs'))
    if 'sort_col' in this_args:
        sort_col = int(this_args.get('sort_col'))
    if 'threshold' in this_args:
        thr = float(this_args.get('threshold'))

    return (n, expression_threshold, sort_col, thr)


@app.route('/top_genes/<column_wildcard>', methods=['GET'])
#@auto.doc()
def gene_index( column_wildcard ):
    corr_type = 'whole-gene'
    usage = 'Example Usage: top_genes/T1T3?n=10&min_exprs=100&sort_col=1'

    expression_threshold=0
    sort_col = 1
    n = int(request.args.get('n'))
    if 'min_exprs' in request.args:
        expression_threshold = int(request.args.get('min_exprs'))
    if 'sort_col' in request.args:
        sort_col = int(request.args.get('sort_col'))
    result_elements = list() 


    some_results = gene_results.filter(regex=re.escape(column_wildcard))[gene_results['Average (quantile normalized)']>expression_threshold]

    selected_column = some_results.columns[sort_col]
    some_results = some_results.sort_values(selected_column, ascending=n>0)
    n = abs(n)
    some_results = some_results.head( n=n ) 
    some_results.columns = [c.replace('_', ' ') for c in some_results.columns] 
    selected_column = some_results.columns[sort_col]

    gene_set_for_search = '['+''.join( ["{'gene':'"+c+"'}," for c in some_results.index])+']'

    scatterize_all_data = ';'.join(some_results.index)
    scatterize_all_link = '<a href="../scatterize_list/genes?list='+scatterize_all_data+'">Scatterize these genes. </a>'
    scatterize_link_notes = 'This will link to a Scatterize page with all these genes along with various behavioral and physiological measures of interst, including AT.'

    enrichr_all_data = ';'.join(some_results.index)
    enrichr_all_label = selected_column.replace(' ','_')+'_top_'+str(n)+'_genes_over_'+str(expression_threshold)+'_reads'
    enrichr_all_link = '<a href="../enrichr_list/'+enrichr_all_data+'?analysis_name='+enrichr_all_label+'">Enrichr these genes. <i class="fa fa-external-link" aria-hidden="true"></i></a>'
    enrichr_link_notes = 'This will link to Enrichr for gene enrichement analyses of the genes listed on this page.' 

    export_all_data = ';'.join(some_results.index)
    export_all_label = selected_column.replace(' ','_')+'_top_'+str(n)+'_genes_over_'+str(expression_threshold)+'_reads'
    export_all_link = '<a href="../export_list/'+export_all_label+'.txt?list='+export_all_data+'">Export this gene list.</a>' 
    export_link_notes = 'This will return a .txt with the genes on this page.' 

    some_results['Gene Name'] = ["<a href=\"/results/"+c+"\">"+c+"</a>" for c in some_results.index ]

    cols = some_results.columns.values
    cols = list(cols[-1:]) + list(cols[:-1])
    some_results = some_results[cols]

    pd.set_option('display.max_colwidth', -1)
    gene_list = some_results.to_html( classes='table table-striped', escape=False, index=False)
    pd.reset_option('display.max_colwidth')


    gene_list_notes = 'The "Gene Name" links to more info on the gene.'
    result_elements.append( {'title': 'Top Gene List ('+column_wildcard+')', 'notes': gene_list_notes, 'content': Markup(gene_list) } )
    result_elements.append( {'title': 'Scatterize feature list', 'notes': scatterize_link_notes, 'content': Markup(scatterize_all_link) } )
    result_elements.append( {'title': 'Enrichr feature list', 'notes': enrichr_link_notes, 'content': Markup(enrichr_all_link) } )
    result_elements.append( {'title': 'Export feature list', 'notes': export_link_notes, 'content': Markup(export_all_link) } )


    this_title = selected_column
    return render_template('top_list.html', **locals())
    

@app.route('/top_genes_from_features/<column_wildcard>', methods=['GET'])
#@auto.doc()
def genes_from_features_index( column_wildcard ):
    corr_type = 'whole-gene from exon'
    usage = 'Example Usage: top_genes_from_features/T1T3?n=10&min_exprs=100&sort_col=1'
    result_elements = list() 

    (n, expression_threshold, sort_col, thr) = get_top_list_args( request.args )  
    some_results = gene_from_features_results.filter(regex=re.escape(column_wildcard))# [gene_from_features_results['Average (quantile normalized)']>expression_threshold]

    some_results = some_results.ix[some_results.filter(regex='_df$').min(axis=1)>1]

    if not n:
        n = some_results.shape[0]

    selected_column = some_results.columns[np.abs(sort_col)]
    some_results = some_results.sort_values(selected_column, ascending=n>0)
    n = abs(n)
    if thr:
        n = min(n, sum(some_results[selected_column]<thr) )

    some_results = some_results.head( n=n ) 
    some_results.columns = [c.replace('_', ' ') for c in some_results.columns] 
    selected_column = some_results.columns[sort_col]

    gene_set_for_search = '['+''.join( ["{'gene':'"+c+"'}," for c in some_results.index])+']'

    scatterize_all_data = ';'.join(some_results.index)
    scatterize_all_link = '<a href="../scatterize_list/genes?list='+scatterize_all_data+'">Scatterize these genes. </a>'
    scatterize_link_notes = 'This will link to a Scatterize page with all these genes along with various behavioral and physiological measures of interst, including AT.'

    enrichr_all_data = ';'.join(some_results.index)
    enrichr_all_label = selected_column.replace(' ','_')+'_top_'+str(n)+'_genes_over_'+str(expression_threshold)+'_reads'
    enrichr_all_link = '<a href="../enrichr_list/'+enrichr_all_data+'?analysis_name='+enrichr_all_label+'">Enrichr these genes. <i class="fa fa-external-link" aria-hidden="true"></i></a>'
    enrichr_link_notes = 'This will link to Enrichr for gene enrichement analyses of the genes listed on this page.' 

    export_all_data = ';'.join(some_results.index)
    export_all_label = selected_column.replace(' ','_')+'_top_'+str(n)+'_genes_over_'+str(expression_threshold)+'_reads'
    export_all_link = '<a href="../export_list/'+export_all_label+'.txt?list='+export_all_data+'">Export this list.</a>' 
    export_link_notes = 'This will return a .txt with the genes on this page.' 

    some_results['Gene Name'] = ["<a href=\"/results/"+c+"\">"+c+"</a>" for c in some_results.index ]

    cols = some_results.columns.values
    cols = list(cols[-1:]) + list(cols[:-1])
    some_results = some_results[cols]

    pd.set_option('display.max_colwidth', -1)
    gene_list = some_results.to_html( classes='table table-striped', escape=False, index=False)
    pd.reset_option('display.max_colwidth')


    gene_list_notes = 'The "Gene Name" links to more info on the gene.'
    result_elements.append( {'title': 'Top Gene List ('+column_wildcard+')', 'notes': gene_list_notes, 'content': Markup(gene_list) } )
    result_elements.append( {'title': 'Scatterize feature list', 'notes': scatterize_link_notes, 'content': Markup(scatterize_all_link) } )
    result_elements.append( {'title': 'Enrichr feature list', 'notes': enrichr_link_notes, 'content': Markup(enrichr_all_link) } )
    result_elements.append( {'title': 'Export feature list', 'notes': export_link_notes, 'content': Markup(export_all_link) } )


    this_title = selected_column
    return render_template('top_list.html', **locals())
    





@app.route('/top_features/<column_wildcard>', methods=['GET'] )
#@auto.doc()
def feature_index( column_wildcard ):
    corr_type = 'feature'
    usage = 'Example Usage: top_features/T1T3?n=10&min_exprs=100&sort_col=1'

    expression_threshold=0
    sort_col = 2
    if 'n' in request.args:
        n = int(request.args.get('n'))
    else: 
        n=results.shape[0]
        print( n )
    if 'min_exprs' in request.args:
        expression_threshold = int(request.args.get('min_exprs'))
    if 'sort_col' in request.args:
        sort_col = int(request.args.get('sort_col'))
    result_elements = list() 

    some_results = results.filter(regex=re.escape(column_wildcard)+'|^name$')[results['Average (quantile normalized)']>expression_threshold]
    selected_column = some_results.columns[np.abs(sort_col)]
    some_results = some_results.sort_values(selected_column, ascending=n>0)

    n = abs(n)
    if 'threshold' in request.args:
        thr = float(request.args.get('threshold'))
        n = min(n, sum(some_results[selected_column]<thr) )

    some_results = some_results.head( n=n ) 

    some_results['Gene Name'] = [c.replace('@',',').split(',')[1] for c in some_results.name ]
    some_results.columns = [c.replace('_', ' ') for c in some_results.columns] 
    selected_column = some_results.columns[sort_col]


    gene_set_for_search = '['+''.join( ["{'gene':'"+c+"'}," for c in some_results['Gene Name']] )+']'
    
    # must be done before adding links... 
    scatterize_all_data = ';'.join(some_results.name)
    scatterize_all_link = '<a href="../scatterize_list/features?list='+scatterize_all_data+'">Scatterize these genes. </a>'
    scatterize_link_notes = 'This will link to a Scatterize page with all these features along with various behavioral and physiological measures of interst, including AT.'

    enrichr_all_data = ';'.join(some_results['Gene Name'])
    enrichr_all_label = selected_column.replace(' ','_')+'_top_'+str(n)+'_features_over_'+str(expression_threshold)+'_reads'
    enrichr_all_link = '<a href="../enrichr_list/'+enrichr_all_data+'?analysis_name='+enrichr_all_label+'">Enrichr these genes. <i class="fa fa-external-link" aria-hidden="true"></i></a>'
    enrichr_link_notes = 'This will link to Enrichr for gene enrichement analyses of the genes listed on this page.' 

    david_gene_id_list = [c.replace(':',',').split(',')[1] for c in some_results.name ]
    david_all_data = ','.join(david_gene_id_list)
    david_all_label =  selected_column.replace(' ','_')+'_top_'+str(n)+'_features_over_'+str(expression_threshold)+'_reads'
    david_all_link = '<a href="../david_list/'+david_all_data+'?analysis_name='+david_all_label+'">David these genes. <i class="fa fa-external-link" aria-hidden="true"></i></a>'
    david_link_notes = 'This will link to David for gene enrichement analyses of the genes listed on this page.' 

    export_all_data = ';'.join(some_results.name)
    export_all_label = selected_column.replace(' ','_')+'_top_'+str(n)+'_features_over_'+str(expression_threshold)+'_reads'
    export_all_link = '<a href="../export_list/'+export_all_label+'.txt?list='+export_all_data+'">Export this list</a>' 
    export_link_notes = 'This will return a .txt with the features on this page.' 

    # convert gene names to links. 
    some_results['Gene Name'] = ["<a href=\"/results/"+c+"\">"+c+"</a>" for c in some_results['Gene Name'] ]
    some_results['name'] = ["<a href=\"/scatterize_feature/"+c+"\">"+c+"</a>" for c in some_results['name'] ]

    cols = some_results.columns.values
    cols = list(cols[-1:]) + list(cols[:-1])
    some_results = some_results[cols]

    pd.set_option('display.max_colwidth', -1)
    feature_list = some_results.to_html( classes='table table-striped', escape=False, index=False)
    pd.reset_option('display.max_colwidth')

    feature_list_notes = 'The "Gene Name" links to more info on the gene and the feature "name" links to a scatterize plot.'
    result_elements.append( {'title': 'Top Feature List ('+column_wildcard+')', 'notes': feature_list_notes, 'content': Markup(feature_list) } )
    result_elements.append( {'title': 'Scatterize feature list', 'notes': scatterize_link_notes, 'content': Markup(scatterize_all_link) } )
    result_elements.append( {'title': 'Enrichr feature list', 'notes': enrichr_link_notes, 'content': Markup(enrichr_all_link) } )
    result_elements.append( {'title': 'David feature list', 'notes': david_link_notes, 'content': Markup(david_all_link) } )
    result_elements.append( {'title': 'Export feature list', 'notes': export_link_notes, 'content': Markup(export_all_link) } )
    this_title = selected_column

    return render_template('top_list.html', **locals())



@app.route('/export_list/<file_name>', methods=['GET'])
#@auto.doc()
def export_list(file_name):
    if 'list' in request.args:
        gene_list = request.args.get('list')
    else:
        gene_list = 'NTRK3;RPS6KA3;APP;CRHR1'
    def generate_list():
        yield '\n'.join( gene_list.split(';') )
    return Response(generate_list(), mimetype='text/csv', headers={"Content-Disposition": "attachment;filename="+file_name} )



@app.route('/scatterize_list/<genes_or_features>', methods=['GET'])
#@auto.doc()
def scatterize_list( genes_or_features ):
    list_for_scatterize = request.args.get('list')
    list_for_scatterize = list_for_scatterize.split(';')

    other_cols = ['Freezing ','Cooing ','Cortisol ','Anxious ','Age ','Not included', 'Relocation']
    this_regex = '^'+'$|^'.join(list_for_scatterize)+'$|'+'|'.join(other_cols)
    this_d = alld.filter(regex=this_regex)

    AT_idx = this_d.columns.get_loc("Anxious Temperament (mean)") + 1

    nus_idx = ''
    nus_idx = nus_idx+str(this_d.columns.get_loc('Age (Time 2)')+1)
    nus_idx = nus_idx+','+str(this_d.columns.get_loc('Not included in Fox et al., 2012')+1)
    nus_idx = nus_idx+','+str(this_d.columns.get_loc('Age when RNA was taken')+1)
    nus_idx = nus_idx+','+str(this_d.columns.get_loc('Relocation Stress')+1)

    url = scatterize_this( this_d )
    return redirect(url+'#x='+str(AT_idx)+'&y=1&n='+nus_idx)
    

def scatterize_this( this_dataframe ):
    my_csv = StringIO()

    this_dataframe.to_csv(my_csv)
    my_csv.seek(0)
    files = {'csvfile': ('for_scatterize.csv', my_csv.read() ) }
    url = 'http://webtasks.keck.waisman.wisc.edu/scatterize/d'
    r = requests.post(url, files=files)
    return( r.url )


@app.route('/scatterize/<gene_of_interest>')
#@auto.doc()
def scatterize( gene_of_interest ):
    other_cols = ['Freezing ','Cooing ','Cortisol ','Anxious ','Age ','Not included', 'Relocation']
    this_regex = '|'.join(other_cols)+'|^'+gene_of_interest+'$|,'+gene_of_interest+'@'
    this_d = alld.filter(regex=this_regex)

    AT_idx = this_d.columns.get_loc("Anxious Temperament (mean)") + 1

    nus_idx = ''
    nus_idx = nus_idx+str(this_d.columns.get_loc('Age (Time 2)')+1)
    nus_idx = nus_idx+','+str(this_d.columns.get_loc('Not included in Fox et al., 2012')+1)
    nus_idx = nus_idx+','+str(this_d.columns.get_loc('Age when RNA was taken')+1)
    nus_idx = nus_idx+','+str(this_d.columns.get_loc('Relocation Stress')+1)

    url = scatterize_this( this_d )
    return redirect( url+'#x='+str(AT_idx)+'&y=1&n='+nus_idx)

@app.route('/scatterize_feature/<feature_of_interest>')
#@auto.doc()
def scatterize_feature( feature_of_interest ):
    other_cols = ['Freezing ','Cooing ','Cortisol ','Anxious ','Age ','Not included', 'Relocation']
    this_regex = '|'.join(other_cols)+'|^'+feature_of_interest+'$'
    this_d = alld.filter(regex=this_regex)

    AT_idx = this_d.columns.get_loc("Anxious Temperament (mean)") + 1

    nus_idx = ''
    nus_idx = nus_idx+str(this_d.columns.get_loc('Age (Time 2)')+1)
    nus_idx = nus_idx+','+str(this_d.columns.get_loc('Not included in Fox et al., 2012')+1)
    nus_idx = nus_idx+','+str(this_d.columns.get_loc('Age when RNA was taken')+1)
    nus_idx = nus_idx+','+str(this_d.columns.get_loc('Relocation Stress')+1)

    url = scatterize_this( this_d )
    return redirect( url+'#x='+str(AT_idx)+'&y=1&n='+nus_idx)


@app.route('/plot_features/<gene_of_interest>.png')
#@auto.doc()
def plot_feature_k(gene_of_interest):
    column_regex=','+gene_of_interest+'@'
    these_results = np.empty([d.filter(regex=column_regex).shape[1],2])
    for i,c in enumerate(d.filter(regex=column_regex).columns):
        x=c
        varlist='^'+x+'$'

        tmp = d.filter(regex=varlist).dropna()
        these_results[i,0] = tmp.mean()
        these_results[i,1] = tmp.std()

    these_results = pd.DataFrame(these_results, columns=['average_exprs', 'std_exprs'])
    these_results.index = d.filter(regex=column_regex).columns

    fig = plt.figure(figsize=(18, 3), facecolor='white')
    ax = fig.add_subplot(111, axisbg='w')

    max_out_colormap=np.ceil(max(these_results['average_exprs']))
    min_out_colormap=np.floor(min(these_results['average_exprs']))
    my_cmap = cm.get_cmap('BuPu') # or any other one such as PiYG
    norm = colors.Normalize(min_out_colormap, max_out_colormap) # the color maps work for [0, 1]

    fig_start = None
    fig_end = None
    for feature in d.filter(regex=column_regex).columns:
        amt = these_results.ix[these_results.index==feature]['average_exprs'].values[0]

        p = re.split(r'[,@\:\=\-\_]', feature)   

        gene_name = p[2]
        feature_type = p[3]
        coords = [int(coord) for coord in p[5:len(p)] ]
        start = 1.0*coords[0]
        end = 1.0*coords[-1]
        middle = start+(end-start)/2.0

        this_color = (0.0,0.0,0.0) # returns an rgba value

        # PLOT JUNCTION
        if feature_type[1:] == 'JXN':
            for c_num,c in enumerate(coords[:-1]):
                ax.annotate("", xy=(c, 20.), xytext=(coords[c_num+1], 20.), textcoords='data',
                    arrowprops=dict(arrowstyle="<->", connectionstyle="angle,angleA=135,angleB=45,rad=0.0",
                    ec=this_color ) )

        # PLOT EXON
        elif feature_type == 'EXON':
            exon = patches.Rectangle((start,-10), end-start, 20, facecolor=this_color, ls='solid', edgecolor='k', linewidth=0.5 )
            ax.add_patch(exon)

        # PLOT INTRON
        elif feature_type == 'ITRN':
            intron = patches.Rectangle((start,-2.5), end-start, 5, color=this_color)
            ax.add_patch(intron)

        if fig_start:
            fig_start = min(start, fig_start)
        else:
            fig_start = start

        if fig_end:
            fig_end = max(end, fig_end)
        else:
            fig_end = end

    plt.title(gene_of_interest+' (length='+str(fig_end-fig_start)+')', fontsize=20, fontstyle='italic', loc='left')

    plt.xlim([fig_start, fig_end])
    x_formatter = ticker.ScalarFormatter(useOffset=False)
    ax.xaxis.set_major_formatter(x_formatter)

    plt.ylim([-100, 100]) #,.1*fig_end-fig_start])
    ax.axes.get_yaxis().set_visible(False)

    cmmapable = cm.ScalarMappable(norm, my_cmap)
    cmmapable.set_array(range(int(min_out_colormap),int(max_out_colormap)))
    cbar = plt.colorbar(cmmapable)
    cbar.set_label('feature expression amount', rotation=270, labelpad=11)

    canvas=FigureCanvas(fig)
    plt.tight_layout()
    png_output = StringIO()
    canvas.print_png(png_output)
    response=make_response(png_output.getvalue())
    response.headers['Content-Type'] = 'image/png'
    return response



@app.route('/plot_features_AT/<gene_of_interest>.png')
def plot_features_AT(gene_of_interest):
    # show the gene
    column_regex=','+gene_of_interest+'@'
    # print column_regex
    these_results = results
    these_results.index = these_results['name']
    these_results = these_results.transpose()
    these_results = these_results.filter(regex=column_regex )

    some_results = these_results.transpose().filter(regex='^Anxious Temperament \(mean\)')
    # print these_results.transpose()['Anxious_Temperament_Time1TimeOD_mean']


    fig = plt.figure(figsize=(18, 3), facecolor='white')
    ax = fig.add_subplot(111, axisbg='w')

    max_out_colormap=4
    my_cmap = cm.get_cmap('RdBu') # or any other one such as PiYG
    norm = colors.Normalize(-1*max_out_colormap, max_out_colormap) # the color maps work for [0, 1]


    fig_start = None
    fig_end = None
    for feature in d.filter(regex=column_regex).columns:
        t_value = some_results.ix[some_results.index==feature].values[0,0]
        p_value = some_results.ix[some_results.index==feature].values[0,1]

        p = re.split(r'[,@\:\=\-\_]', feature)   

        gene_name = p[2]
        feature_type = p[3]
        coords = [int(coord) for coord in p[5:len(p)] ]
        start = 1.0*coords[0]
        end = 1.0*coords[-1]
        middle = start+(end-start)/2.0

    #     print gene_name, feature_type, t_value, start, end, coords

        rank_t = max( min(t_value, max_out_colormap), -1*max_out_colormap)
        this_color = my_cmap(norm(rank_t*-1)) # returns an rgba value
    #     print this_color

        # PLOT JUNCTION
        if feature_type[1:] == 'JXN':
            for c_num,c in enumerate(coords[:-1]):
                ax.annotate("", xy=(c, 20.), xytext=(coords[c_num+1], 20.), textcoords='data',
                    arrowprops=dict(arrowstyle="<->", connectionstyle="angle,angleA=135,angleB=45,rad=0.0",
                    ec=this_color ) )
            if p_value<.05:
                ax.annotate(round(t_value,2),xy=(start, 40.), xytext=(middle, 40.) )

        # PLOT EXON
        elif feature_type == 'EXON':
            exon = patches.Rectangle((start,-10), end-start, 20, color=this_color )
            ax.add_patch(exon)
            if p_value<.05:
                ax.annotate(round(t_value,2),
                    xy=(start, 12.), xytext=(middle, 12.), fontsize=20 )

        # PLOT INTRON
        elif feature_type == 'ITRN':
            intron = patches.Rectangle((start,-2.5), end-start, 5, color=this_color)
            ax.add_patch(intron)
            if p_value<.05:
                ax.annotate(round(t_value,2),
                    xy=(start, -20.0), xytext=(middle, -20.0))

        if fig_start:
            fig_start = min(start, fig_start)
        else:
            fig_start = start

        if fig_end:
            fig_end = max(end, fig_end)
        else:
            fig_end = end



    plt.title(gene_of_interest+' (length='+str(fig_end-fig_start)+')', fontsize=20, fontstyle='italic', loc='left')

    plt.xlim([fig_start, fig_end])
    x_formatter = ticker.ScalarFormatter(useOffset=False)
    ax.xaxis.set_major_formatter(x_formatter)

    plt.ylim([-100, 100]) #,.1*fig_end-fig_start])
    ax.axes.get_yaxis().set_visible(False)

    cmmapable = cm.ScalarMappable(norm, my_cmap)
    cmmapable.set_array(range(-1*max_out_colormap,max_out_colormap))
    cbar = plt.colorbar(cmmapable)
    cbar.set_label('feature t-value', rotation=270, labelpad=11)


    # plt.savefig(gene_of_interest+'.png')

    canvas=FigureCanvas(fig)
    plt.tight_layout()
    png_output = StringIO()
    canvas.print_png(png_output)
    response=make_response(png_output.getvalue())
    response.headers['Content-Type'] = 'image/png'
    return response


@app.route('/plot_feature_scatters/<feature_of_interest>.png')
def plot_feature_scatters( feature_of_interest, drop_x_greaterthan=None, drop_x_lessthan=None, outfile=None, rank=None ):
    x = feature_of_interest
    if drop_x_greaterthan:
        this_d=alld[ alld[x]<drop_x_greaterthan]
    elif drop_x_lessthan:
        this_d=alld[ alld[x]>drop_x_lessthan]
    else:
        this_d = alld
        
    if rank:
        this_d = this_d.rank()
        
    f, axarr = plt.subplots(3,3)
#     f.set_size_inches(18.5,10.5)
    f.set_size_inches(11,11)
    f.set_facecolor('white')

    y='Freezing_duration_Time1'
    ax = sns.regplot(x,y,this_d.filter([x,y]).dropna(), ax=axarr[0,0])
    ax.set_xlabel(x.replace(',','\n').replace(':COORDS=','\n'))
    y='Cooing_frequency_Time1'
    ax = sns.regplot(x,y,this_d.filter([x,y]).dropna(), ax=axarr[0,1])
    ax.set_xlabel(x.replace(',','\n').replace(':COORDS=','\n'))
    y='Cortisol_levels_Time1'
    ax = sns.regplot(x,y,this_d.filter([x,y]).dropna(), ax=axarr[0,2])
    ax.set_xlabel(x.replace(',','\n').replace(':COORDS=','\n'))

    
    y='Freezing_duration_Time2'
    ax = sns.regplot(x,y,this_d.filter([x,y]).dropna(), ax=axarr[1,0])
    ax.set_xlabel(x.replace(',','\n').replace(':COORDS=','\n'))
    y='Cooing_frequency_Time2'
    ax = sns.regplot(x,y,this_d.filter([x,y]).dropna(), ax=axarr[1,1])
    ax.set_xlabel(x.replace(',','\n').replace(':COORDS=','\n'))
    y='Cortisol_levels_Time1'
    ax = sns.regplot(x,y,this_d.filter([x,y]).dropna(), ax=axarr[1,2])
    ax.set_xlabel(x.replace(',','\n').replace(':COORDS=','\n'))

    y='Anxious_Temperament_Time1Time2_mean'
    ax = sns.regplot(x,y,this_d.filter([x,y]).dropna(), ax=axarr[2,0])
    ax.set_xlabel(x.replace(',','\n').replace(':COORDS=','\n'))
    y='AnxTemp_mean_T1T3'
    ax = sns.regplot(x,y,this_d.filter([x,y]).dropna(), ax=axarr[2,1])
    ax.set_xlabel(x.replace(',','\n').replace(':COORDS=','\n'))
    y='Anxious_Temperament_Time1ToD_clust_pos_23'
    ax = sns.regplot(x,y,this_d.filter([x,y]).dropna(), ax=axarr[2,2])
    ax.set_xlabel(x.replace(',','\n').replace(':COORDS=','\n'))
    
    axarr[2,2].set_ylim([np.min(this_d[y]),np.max(this_d[y])])

    canvas=FigureCanvas(f)
    plt.tight_layout()
    png_output = StringIO()
    canvas.print_png(png_output)
    response=make_response(png_output.getvalue())
    response.headers['Content-Type'] = 'image/png'
    return response


@app.route('/feature_list_<gene_of_interest>')
def feature_list(gene_of_interest):
    column_regex=','+gene_of_interest+'@'
    return '%s' % d.filter(regex=column_regex).columns.values

@app.route('/plot_genePET_scatters/<gene_of_interest>.png')
def plot_genePET_scatters( gene_of_interest, drop_x_greaterthan=None, drop_x_lessthan=None, outfile=None, rank=None ):
    x = gene_of_interest
    if drop_x_greaterthan:
        this_d=gened[ gened[x]<drop_x_greaterthan]
    elif drop_x_lessthan:
        this_d=gened[ gened[x]>drop_x_lessthan]
    else:
        this_d = gened
        
    if rank:
        this_d = this_d.rank()
        
    f, axarr = plt.subplots(3,3)
#     f.set_size_inches(18.5,10.5)
    f.set_size_inches(11,11)
    f.set_facecolor('white')

    y='Anxious_Temperament_Time1ToD_clust_pos_23'
    sns.regplot(x,y,this_d.filter([x,y]).dropna(), ax=axarr[0,0])
    y='Anxious_Temperament_Time1ToD_clust_pos_15'
    sns.regplot(x,y,this_d.filter([x,y]).dropna(), ax=axarr[0,1])
    y='Anxious_Temperament_Time1ToD_clust_pos_16'
    sns.regplot(x,y,this_d.filter([x,y]).dropna(), ax=axarr[0,2])

    
    y='Anxious_Temperament_Time1ToD_clust_pos_8'
    sns.regplot(x,y,this_d.filter([x,y]).dropna(), ax=axarr[1,0])
    y='Anxious_Temperament_Time1ToD_clust_pos_9'
    sns.regplot(x,y,this_d.filter([x,y]).dropna(), ax=axarr[1,1])
    y='Anxious_Temperament_Time1ToD_clust_pos_7'
    sns.regplot(x,y,this_d.filter([x,y]).dropna(), ax=axarr[1,2])

    y='Anxious_Temperament_Time1ToD_clust_pos_1'
    sns.regplot(x,y,this_d.filter([x,y]).dropna(), ax=axarr[2,0])
    y='Anxious_Temperament_Time1ToD_clust_pos_2'
    sns.regplot(x,y,this_d.filter([x,y]).dropna(), ax=axarr[2,1])
    y='Anxious_Temperament_Time1ToD_clust_pos_3'
    sns.regplot(x,y,this_d.filter([x,y]).dropna(), ax=axarr[2,2])
    
    axarr[2,2].set_ylim([np.min(this_d[y]),np.max(this_d[y])])
    
    canvas=FigureCanvas(f)
    plt.tight_layout()
    png_output = StringIO()
    canvas.print_png(png_output)
    response=make_response(png_output.getvalue())
    response.headers['Content-Type'] = 'image/png'
    return response


@app.route('/plot_gene_scatters/<gene_of_interest>.png')
#@auto.doc()
def plot_gene_scatters( gene_of_interest, drop_x_greaterthan=None, drop_x_lessthan=None, outfile=None, rank=None ):
    x = gene_of_interest
    if drop_x_greaterthan:
        this_d=gened[ gened[x]<drop_x_greaterthan]
    elif drop_x_lessthan:
        this_d=gened[ gened[x]>drop_x_lessthan]
    else:
        this_d = gened
        
    if rank:
        this_d = this_d.rank()
        
    f, axarr = plt.subplots(3,3)
#     f.set_size_inches(18.5,10.5)
    f.set_size_inches(11,11)
    f.set_facecolor('white')

    y='Freezing_duration_Time1'
    sns.regplot(x,y,this_d.filter([x,y]).dropna(), ax=axarr[0,0])
    y='Cooing_frequency_Time1'
    sns.regplot(x,y,this_d.filter([x,y]).dropna(), ax=axarr[0,1])
    y='Cortisol_levels_Time1'
    sns.regplot(x,y,this_d.filter([x,y]).dropna(), ax=axarr[0,2])

    
    y='Freezing_duration_Time2'
    sns.regplot(x,y,this_d.filter([x,y]).dropna(), ax=axarr[1,0])
    y='Cooing_frequency_Time2'
    sns.regplot(x,y,this_d.filter([x,y]).dropna(), ax=axarr[1,1])
    y='Cortisol_levels_Time1'
    sns.regplot(x,y,this_d.filter([x,y]).dropna(), ax=axarr[1,2])

    y='Anxious_Temperament_Time1Time2_mean'
    sns.regplot(x,y,this_d.filter([x,y]).dropna(), ax=axarr[2,0])
    y='AnxTemp_mean_T1T3'
    sns.regplot(x,y,this_d.filter([x,y]).dropna(), ax=axarr[2,1])
    y='ATPfcRCONN_robust'
    sns.regplot(x,y,this_d.filter([x,y]).dropna(), ax=axarr[2,2])

    canvas=FigureCanvas(f)
    plt.tight_layout()
    png_output = StringIO()
    canvas.print_png(png_output)
    response=make_response(png_output.getvalue())
    response.headers['Content-Type'] = 'image/png'
    return response


def format_t_p_table(df):
    p = df.filter(regex='_p$').transpose()
    p.index = [c[:-2] for c in p.index]
    t = df.filter(regex='_t$').transpose()
    t.index = [c[:-2] for c in t.index]
    df = t.merge(p, left_index=True, right_index=True )
    df.columns = ['t-value', 'p-value']

    df.index = [c.replace('_',' ') for c in df.index ]

    significant = lambda x: '<span class="significant_text">%1.6f</span>' % x if x<0.05 else '%1.6f'%x
    df_html = df.to_html(float_format=lambda x:'%1.6f'%x, formatters={'p-value': significant},  classes='table table-striped', escape=False)

    return df_html

def format_table(df, list_of_col_regex, list_of_titles):
    table_df = pd.DataFrame()
    for col in list_of_col_regex:
        this = df.filter(regex=col).transpose()
        this.index = [re.sub(col,'',c) for c in this.index]
        table_df = table_df.merge(this, left_index=True, right_index=True, how='outer')
    table_df.columns = list_of_titles
    significant = lambda x: '<span class="significant_text">%1.6f</span>' % x if x<0.05 else '%1.6f'%x

    df_html = table_df.to_html(float_format=lambda x:'%1.6f'%x, formatters={'p-value': significant},  classes='table table-striped', escape=False)
    return df_html




@app.route('/results', )
def redirect_to_print_results(methods=['GET']):
    print( request.args )
    if 'gene_search' in request.args:
        gene_name = request.args.get('gene_search')
    else:
        gene_name = 'CRH'
    return redirect('results/'+gene_name)

@app.route('/results/')
def results_fail():
    return redirect('error')


@app.route('/results/<gene_of_interest>')
#@auto.doc()
def print_results(gene_of_interest, primary_variable='Anxious Temperament (mean)'):
    result_elements = list()
    if gene_of_interest not in gene_results.index:
        return redirect('error?error_text='+gene_of_interest)


    column_regex=','+gene_of_interest+'@'
    feature_results = results
    feature_results.index = feature_results['name']
    feature_results = feature_results.transpose()
    feature_results = feature_results.filter(regex=column_regex )

    some_results = feature_results.transpose().filter(regex='^'+re.escape(primary_variable))
    some_results.columns = [c.replace(primary_variable, '') for c in some_results.columns]
    some_results['Feature Name'] = ["<a href=\"/scatterize_feature/"+c+"\">"+c+"</a>" for c in some_results.index ]
    some_results.columns = [c.replace('_', ' ') for c in some_results.columns]
    # some_results.columns = ['t-value', 'p-value', 'Feature Name']
    some_results.columns = [c.replace('p', 'p-value') for c in some_results.columns]
    some_results.columns = [c.replace('^t', 't-value') for c in some_results.columns]
    some_results.columns = [c.replace(' ', '') for c in some_results.columns]
    some_results.index = [c.replace('@', ' ') for c in some_results.index]

    cols = some_results.columns.values
    cols = list(cols[-1:]) + list(cols[:-1])
    some_results = some_results[cols]


    # # should replace with a call to format_table
    significant = lambda x: '<span class="significant_text">%1.6f</span>' % x if x<0.05 else '%1.6f'%x
    some_results['p-value'] = [significant(p) for p in some_results['p-value'] ]

    pd.set_option('display.max_colwidth', -1)
    feature_result_content = some_results.to_html(classes='table table-striped', escape=False, index=False)
    pd.reset_option('display.max_colwidth')


    this_title = gene_of_interest

    # gene_feature_img = markdown.markdown("[![gene_model](../plot_features/"+gene_of_interest+".png)](../plot_features/"+gene_of_interest+".png)")
    # gene_feature_AT_img = markdown.markdown("[![gene_model](../plot_features_AT/"+gene_of_interest+".png)](../plot_features_AT/"+gene_of_interest+".png)")
    # gene_scatter_img = markdown.markdown("[scatter plot for gene vs. AT vars](../plot_gene_scatters/"+gene_of_interest+".png)")
    # gene_scatterPET_img = markdown.markdown("[scatter plots for gene vs. PET vars](../plot_genePET_scatters/"+gene_of_interest+".png)")

    gene_result_content = format_t_p_table(gene_results[gene_results.index==gene_of_interest])
    # gene_result_content.index = [c.replace('_',' ') for c in gene_result_content.index ]
    # gene_result_content = gene_result_content.to_html( classes='table table-striped')


    # WOULD NEED TO BE DONE FOR EACH FEATURE
    # feature_scatter_img = markdown.markdown("[scatter plot for feature vs. AT vars](../plot_feature_scatters/"+gene_of_interest+".png)")

    # content = Markup(gene_result_content+text+img+feature_result_content)

    rhesus2human_gene_result_mean =  rhesus2human_gene_results[
        rhesus2human_gene_results.index==gene_of_interest][[
        'Average (quantile normalized)', 
        'Standard deviation', 'Observed in __ subjects'
        ]].to_html( classes='table table-striped' )
    rhesus2human_gene_result_mean_notes = Markup(markdown.markdown('Reads for gene-level data after aligning to the human genome.'))

    rhesus2human_gene_result_content = format_t_p_table(rhesus2human_gene_results[rhesus2human_gene_results.index==gene_of_interest])
    rhesus2human_gene_result_notes = Markup(markdown.markdown('Gene-level associations after aligning to the human genome.'))


    try:
        gene_from_features_result_content = format_table(gene_from_features_results[gene_from_features_results.index==gene_of_interest], ['_R2$', '_df$', '_F_p$'], ['R^2', 'df', 'p-value'])
        gene_from_features_result_notes = Markup(markdown.markdown('Gene-level associations when using multiple regression with all Exons as predictors...'))
    except:
        gene_from_features_result_content = ''
        gene_from_features_result_notes = ''


    try:
        gene_from_features_pet_result_content = format_table(gene_from_features_pet_results[gene_from_features_pet_results.index==gene_of_interest], ['_R2$', '_df$', '_F_p$'], ['R^2', 'df', 'p-value'])
        gene_from_features_pet_result_notes = Markup(markdown.markdown('Gene-level associations when using multiple regression with all Exons as predictors...'))
    except:
        gene_from_features_pet_result_content = ''
        gene_from_features_pet_result_notes = ''



    gene_result_mean =  gene_results[gene_results.index==gene_of_interest][[ 'Average (quantile normalized)', 'Standard deviation', 'Observed in __ subjects']
].to_html( classes='table table-striped' )

    gene_PET_result_content = format_t_p_table(gene_pet_results[gene_pet_results.index==gene_of_interest].filter(regex='clust'))

    gene_PET_result_notes = Markup(markdown.markdown('FDG-PET clusters where metabolism was correlated with AT (p<.005). Clusters are numbered from largest to smallest; positive & negative effects are listed seperately. [Click to see cluster-map.](../static/neuroviewer/index_line.html)'))

    result_elements.append( {'title': 'Expression Level', 'content': Markup(gene_result_mean) } )
    #result_elements.append( {'title': 'Gene Model', 'content': Markup(gene_feature_img) } )
    result_elements.append( {'title': 'Gene Result', 'content': Markup(gene_result_content) } )
#    result_elements.append( {'title': 'Gene Model and AT', 'content': Markup(gene_feature_AT_img) } )
    result_elements.append( {'title': 'Gene Model and AT', 'content': Markup('<div id="gene_model_js"> </div>') } )
    result_elements.append( {'title': 'Feature Result', 'notes':'"Feature Name" link goes to scatterize.', 'content': Markup(feature_result_content) } )

    result_elements.append( {'title': 'Gene results from features', 'notes': gene_from_features_result_notes, 'content': Markup(gene_from_features_result_content) } )
    result_elements.append( {'title': 'Expression level (mapped to Human)', 'notes': rhesus2human_gene_result_mean_notes, 'content': Markup(rhesus2human_gene_result_mean) } )
    result_elements.append( {'title': 'Gene Results (mapped to Human)', 'notes': rhesus2human_gene_result_notes, 'content': Markup(rhesus2human_gene_result_content) } )
    result_elements.append( {'title': 'Gene PET results', 'notes': gene_PET_result_notes, 'content': Markup(gene_PET_result_content) } )
    result_elements.append( {'title': 'Gene PET results from features', 'notes': gene_from_features_pet_result_notes, 'content': Markup(gene_from_features_pet_result_content) } )
    # result_elements.append( {'title': 'Gene Level Scatters', 'content': Markup(gene_scatter_img) } )
    # result_elements.append( {'title': 'Gene Level PET Scatters', 'content': Markup(gene_scatterPET_img) } )

    return render_template('results.html', **locals())


def init():
    rnaseq_file = 'static/data/feature_quantification/rhesus_features_and_intergenes/RHESUS_QUANTILE_FEATURES.scrs'
    column_names = ['name','1', '11', '12', '13', '14', '15', '16', '17', '18', '19', '2', '21', '22', '23', '24', '25', '26', '27', '28', '29', '3', '30', '31', '32', '33', '34', '35', '36', '37', '38', '39', '4', '40', '41', '42', '43',     '44', '45', '46', '47', '48', '5', '6', '7', '8', '9', 'type', 'n', 'mean?']
    # read data -- surprisingly hard, make sure to skip the first line... 
    d = pd.read_table(rnaseq_file, delim_whitespace=True, skiprows=1, header=None, names=column_names )
    d.index = d.name
    d = d.ix[:,1:47].copy()
    d = d.transpose()
    d.index = [int(idx) for idx in d.index]

    results = pd.read_csv('static/data/RNAseq_quants_by_feature_ols_quantile_quantification_covAge2AodNorsStress.csv')


    # read gene-level-data
    rnaseq_dir = 'static/data/gene_quantifications/quantile/'
    gene_data = pd.DataFrame()
    files = glob.glob(rnaseq_dir+'/*')
    for f in files:
        cur_id = re.search('static/data/gene_quantifications/quantile/Rh(\d+)\.gene\.quantile', f).group(1)
        cur_file = pd.read_csv(f, index_col=0, header=0, names=['id', int(cur_id)], dtype={int(cur_id): np.float64} , sep=' ')
        gene_data = gene_data.merge(cur_file, left_index=True, right_index=True, how='outer')
    gene_data.index = ['_'.join(s.split(',')[1:]) for s in gene_data.index ]
    gene_data=gene_data.transpose()

    gene_results = pd.read_csv('static/data/RNAseq_quants_by_gene_ols_gene_quantification_covAge2AodNorsStress.csv')
    gene_results.index = gene_results['Unnamed: 0']
    gene_results.index.name = None
    gene_pet_results = pd.read_csv('static/data/RNAseq_quants_by_gene_ols_gene_quantification_clusters_PETT1ToD_ATT1ToD_covAgeT2AgeToDNorsStress.csv')
    gene_pet_results.index = gene_pet_results['Unnamed: 0']
    # gene_results = gene_results.merge( gene_pet_results.filter(regex='clust'), left_index=True, right_index=True )


    gene_from_features_results = pd.read_csv('static/data/RNAseq_FeaturesCombined_quants_by_features_ols_feature_quantification_covAge2AodNorsStress.csv')
    gene_from_features_results.index = gene_from_features_results['Unnamed: 0']

    # print gene_results.columns
    #gene_results.columns[0] = 'Gene Name'
    #gene_results.index = gene_results['Gene Name']

    gene_from_features_pet_results = pd.read_csv('static/data/PET_RNAseq_FeaturesCombined_quants_by_features_ols_feature_quantification_covAge2AodNorsStress.csv')
    gene_from_features_pet_results.index = gene_from_features_pet_results['Unnamed: 0']

    rhesus2human_gene_results = pd.read_csv('static/data/RNAseq_rhesus2human_quants_by_gene_ols_gene_quantification_covAge2AodNorsStress.csv')
    rhesus2human_gene_results.index = rhesus2human_gene_results['Unnamed: 0']
    rhesus2human_gene_results.index.name = None


    phen = pd.read_csv('static/data/WisconsinPhenotypes_Fall2014.csv')
    phen = phen.replace('',np.nan)
    phen.index = phen['USC ID ']
    
    conn_file = 'static/data/connectivity_vals.csv'
    conn = pd.read_csv(conn_file, index_col=0)
    
    setup_file = 'static/data/setup_for_condor.csv'
    setup = pd.read_csv(setup_file, index_col=1)
    ToD = setup[['AT_ToD', 'Cooing_ToD', 'Cortisol_ToD', 'Freezing_ToD', 'AT_mean_T1ToD']]
    
    extracted_PET_T1ToD_data_file = 'static/data/AT_mean_T1ToD_PET_T1ToD_0_t_0025_cluster_over2mm_values_for_pandas.csv'
    extracted_PET_T1ToD_data = pd.read_csv(extracted_PET_T1ToD_data_file)
    extracted_PET_T1ToD_data['MRI_ID_T1'] = [int(c.split('_')[4]) for c in extracted_PET_T1ToD_data.PET_T1ToD]
    extracted_PET_T1ToD_data['MRI_ID_ToD'] = [int(c.split('_')[5]) for c in extracted_PET_T1ToD_data.PET_T1ToD]
    
    this_setup = setup.filter(regex='^MRI_ID').copy()
    this_setup.loc[:,'RNA ID'] = setup.index 
    
    extracted_PET_T1ToD_data = this_setup.merge(extracted_PET_T1ToD_data, on=['MRI_ID_T1', 'MRI_ID_ToD'], how='left' )
    extracted_PET_T1ToD_data.index = extracted_PET_T1ToD_data['RNA ID']
    extracted_PET_T1ToD_data = extracted_PET_T1ToD_data.filter(regex='^Anxious_Temperament_Time1ToD_clust')
    

    alld = d.merge(phen, left_index=True, right_index=True, how='outer')
    alld = gene_data.merge(alld, left_index=True, right_index=True, how='outer')
    # alld = alld.merge(conn, left_on='Subject', right_index=True)
    alld = alld.merge(ToD, left_on='USC ID ', right_index=True)
    alld = alld.merge(extracted_PET_T1ToD_data, left_on='USC ID ', right_index=True)
    
    gened = gene_data.merge(phen, left_index=True, right_index=True, how='outer')
    # gened = gened.merge(conn, left_on='Subject', right_index=True)
    gened = gened.merge(ToD, left_on='USC ID ', right_index=True)
    gened = gened.merge(extracted_PET_T1ToD_data, left_on='USC ID ', right_index=True)


    # THIS IS HOW I SHOULD RENAME
    df_list = [alld, gened,results, gene_results,gene_pet_results, rhesus2human_gene_results,gene_from_features_pet_results,gene_from_features_results ]
    replace_names_dict = {
        '^mean': 'Average (quantile normalized)',
        'raw_mean': 'Average (raw reads)',
        'std': 'Standard deviation',
        'Observed_in_n_subjects': 'Observed in __ subjects',
        'Anxious_Temperament_Time1ToD_mean':'Anxious Temperament (mean)',
        'Anxious_Temperament_Time1TimeOD_mean':'Anxious Temperament (mean)',
        'Anxious_Temperament_Time1':'Anxious Temperament (Time 1)',
        'Anxious_Temperament_Time2':'Anxious Temperament (Time 2)', 
        'Freezing_duration_Time1':'Freezing duration (Time 1)',
        'Freezing_duration_Time2':'Freezing duration (Time 2)',
        'Cooing_frequency_Time1':'Cooing frequency (Time 1)',
        'Cooing_frequency_Time2':'Cooing frequency (Time 2)',
        'Cortisol_levels_Time1':'Cortisol levels (Time 1)',
        'Cortisol_levels_Time2':'Cortisol levels (Time 2)',
        'age_ToD': 'Age when RNA was taken',
        'age_T1': 'Age (Time 1)',
        'age_T2': 'Age (Time 2)',
        'isNORS': 'Not included in Fox et al., 2012',
        'AT_mean_T1ToD': 'Anxious Temperament (mean)',
        'Freezing_mean_T1ToD': 'Freezind duration (mean)',
        'Cooing_mean_T1ToD': 'Cooing frequency (mean)',
        'Cortisol_mean_T1ToD': 'Cortisol Levels (mean)',
        'stress_Group': 'Relocation Stress'
        }


    #[ 'Average (quantile normalized)', 'Average (raw reads)', 'Standard deviation', 'Observed in __ subjects']

    column_keys_to_drop = ['ATPfcRCONN', 'ATUncinateFA', 'Time1Time2_mean']
    for this_df in df_list:
        for col_to_drop in column_keys_to_drop:
            this_df.drop(axis='columns', labels=list(this_df.filter(regex=col_to_drop).columns), inplace=True ) 
        for key in replace_names_dict:
            this_df.columns = this_df.columns.str.replace(key, replace_names_dict[key])
            # this_df.index = this_df.index.str.replace(key, replace_names_dict[key])

    # app.config['SITEMAP_INCLUDE_RULES_WITHOUT_PARAMS'] = True
    # ext.init_app(app)

    return d, alld, gened,results, gene_results,gene_pet_results, rhesus2human_gene_results,gene_from_features_pet_results,gene_from_features_results


d, alld, gened,results, gene_results,gene_pet_results, rhesus2human_gene_results,gene_from_features_pet_results,gene_from_features_results = init()
app.wsgi_app = ProxyFix(app.wsgi_app)

# if __name__ == 'mini_flask_RNAseq_AT' or __name__ == '__main__':
    # app.run(host='127.0.0.1',port=5001)

if __name__ == '__main__':
    app.run()
