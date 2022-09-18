import numpy as np
from scipy import sparse
from scipy.sparse import csr_matrix
import pandas as pd
import time
import random
import rpy2.robjects as robjects
from rpy2.robjects.packages import importr
from rpy2.robjects import r
import os
import getpass
import traceback
import math
import cairosvg
import markov_clustering as mc
import itertools
import IPython
import py4cytoscape as p4c
import networkx as nx
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from pathlib import Path
from io import StringIO
import shutil

def log(i):
	import math
	return -(math.log10(i))

######################################
#####  Bubble plot and Bar graph  #####
######################################
def plt_bubble_bar_svg(pathway_result,outputbubble,outputbar,num,database):
	import matplotlib
	from io import StringIO
	import re
	import pandas as pd
	matplotlib.use('agg')
	import matplotlib.pyplot as plt
	from mpl_toolkits.axes_grid1.inset_locator import inset_axes
	import random
	import os
	import math
	num=int(num)
	data = pd.read_csv(pathway_result,sep='\t')
	data_temp = data.sort_values(by = ['FDR'],axis = 0,ascending = True)[:num]
	data_sort = data_temp.sort_values(by = ['PS'],axis = 0,ascending = True)
	score = list(data_sort['PS'])
	fdr = list(data_sort['FDR'])
	fdr_score=list(map(log,fdr))
	x=list(data_sort['PS'])
	if database=="All":
		y=list(data_sort['pathway_name']+"_"+data_sort['Source'])
	else:
		y=list(data_sort['pathway_name'])

	c=list(fdr_score)
	s = list(data_sort['node_number'])
	fig_len=3.6
	legend_len=2.4
	pathway_len_70= 6
	pathway_len_list=list(map(len,y))
	pathway_len=round((max(pathway_len_list)/70)*6)
	fig_len_total=fig_len+legend_len+pathway_len
	fig_height=7
	if num >20:
		fig_height=num/20*7

	fig, ax = plt.subplots(figsize=(fig_len_total,fig_height))
	axins = inset_axes(ax,width="6%",height="30%",loc='lower left',bbox_to_anchor=(1.1, 0.65, 1, 1),bbox_transform=ax.transAxes)# width = 5% of parent_bbox width,# height : 50%
	import matplotlib.colors
	cmap = matplotlib.colors.LinearSegmentedColormap.from_list("", ["#07F800","#FF0000"],N=512)
	scatter = ax.scatter(x, y, c=c, s=s,cmap=cmap)
	ax.tick_params(axis='both', labelsize=10)
	ax.set_xlabel('PS', fontsize=18)
	#ax.set_title('Pathway Enrichment',fontsize=18)
	legend2 = ax.legend(*scatter.legend_elements(prop='sizes', num = 6),fontsize=10,loc='lower left',title='Gene Number',title_fontsize=12,bbox_to_anchor=(1.01, 0.1),bbox_transform=ax.transAxes,labelspacing=1,frameon=False) 
	clb=fig.colorbar(scatter,shrink=0.3,aspect=20,cax=axins) 
	clb.ax.set_title('-log10(FDR)',loc="center",fontsize=12)
	clb.ax.tick_params(labelsize=10)
	plt.gcf().subplots_adjust(left=(pathway_len/fig_len_total), bottom=0.1,top=0.9,right=((fig_len+pathway_len)/fig_len_total))
	buf = StringIO()
	plt.savefig(fname=buf,format="svg")
	all=buf.getvalue()
	buf.seek(0,0) 
	buf.readline()
	buf.readline()
	buf.readline()
	buf.readline()
	svg=buf.readline() 
	svg_temp=re.sub(r'height=".*?" ',"",svg) 
	svg_out=re.sub(r'width=".*?" ',"",svg_temp) 
	p=re.split('<svg.*?>',all) 
	out=StringIO()
	out.write(p[0])
	out.write(svg_out)
	out.write(p[1]) 
	fd = open (outputbubble, 'w')
	fd.write (out.getvalue())
	fd.close()
	buf.close()
	out.close()
	#############bar plot
	data_temp = data[:num]
	c=list(reversed(list(map(log,data_temp['FDR']))))
	if database=="All":
		y=list(reversed(data_temp['pathway_name']+"_"+data_temp['Source']))
	else:
		y=list(reversed(data_temp['pathway_name']))

	fig, ax = plt.subplots(figsize=(fig_len_total,fig_height))
	fdr_axis_start = math.floor(min(c)*10)/10
	fdr_axis_end = math.ceil(max(c)*10)/10
	plt.xlim(fdr_axis_start,fdr_axis_end)
	ax.barh(y,c,color='#E0683D')
	ax.tick_params(axis='both', labelsize=10) 
	ax.set_xlabel('-log10(FDR)', fontsize=18)#fontstyle='italic'
	plt.gcf().subplots_adjust(left=(pathway_len/fig_len_total), bottom=0.1,top=0.9,right=((fig_len+pathway_len+legend_len/2)/fig_len_total))
	buf = StringIO()
	plt.savefig(fname=buf,format="svg")
	all=buf.getvalue()
	plt.close()
	buf.seek(0,0)
	buf.readline()
	buf.readline()
	buf.readline()
	buf.readline()
	svg=buf.readline()
	svg_temp=re.sub(r'height=".*?" ',"",svg)
	svg_out=re.sub(r'width=".*?" ',"",svg_temp)
	p=re.split('<svg.*?>',all)
	out=StringIO()
	out.write(p[0])
	out.write(svg_out)
	out.write(p[1]) #写入新的svg
	fd = open (outputbar, 'w')
	fd.write (out.getvalue())
	fd.close()
	buf.close()
	out.close()
	plt.close()

############################
#####  Enrichment map  #####
############################
def cytoscape(index,sig_method,sig_value,database,dir):
	import pandas as pd
	import numpy as np
	import markov_clustering as mc
	import itertools
	import networkx as nx
	import py4cytoscape as p4c
	from io import StringIO
	import time
	import cairosvg
	index=float(index)
	#index=0.3
	#sig_value=0.01
	#sig_method="p_value"
	index_use=index
	outdir= dir+"/output"
	outfile= outdir+"/network_temp_"+str(index)
	outtemp=outfile+".svg"
	outnetwork= outdir+"/network_"+str(index)+".svg"
	outpathwayfile=outdir+"/top20pathwayname_"+str(index)+".svg"
	topfile=outdir+'/top20pathway_'+str(index)+'.txt'
	infofile=outdir+'/cluster_info_'+str(index)+'.txt'

	infile= dir+"/output/pathway_result.txt"
	input = pd.read_csv(infile,dtype={'p_value': float},sep="\t")
	index_file=dir+'/data/all_pathway_index_result_matrix.txt'
	sim = pd.read_csv(index_file,sep="\t",index_col=0)
	if sig_method=="FDR":
		used= input.loc[input['FDR']< float(sig_value)]
	else:
		used= input.loc[input['p_value']< float(sig_value)]

	if used.shape[0]==0:
		return

	name=(used['pathway_name']+ "_"+used['Source']).tolist()
	column=sim[name] 
	row=column.loc[name] 
	matrix=row.values
	new_matrix=matrix.copy()
	new_matrix[new_matrix<index]=0
	new_matrix[new_matrix>=index]=1 #
	np.fill_diagonal(new_matrix,0) 
	result = mc.run_mcl(new_matrix) 
	clusters = mc.get_clusters(result) 
	min_p=list()
	for clu in clusters:
		min_p.append(min(used['p_value'][list(clu)])) 

	min_index=np.array(min_p).argsort()[:20] 
	pathway=list() 
	min_pathway=list() 
	node_color=list()
	colors=['#E41A1C','#377EB8','#4DAF4A','#984EA3','#FF7F00','#FFFF33','#A65628','#F781BF','#BD9E39','#8DD3C7','#FFFFB3','#BEBADA',	'#FB8072','#80B1D3','#FDB462','#B3DE69','#FCCDE5','#FDD0A2','#BC80BD','#CCEBC5']
	i=0
	node_pvalue=list()
	cluster_id=list() #1,2,3,...20
	for index in min_index:
		#print(index)
		clu_p=used['p_value'][list(clusters[index])].tolist() 
		min_clu_index=np.array(clu_p).argsort()[:10] 
		path=np.array(clusters[index])[list(min_clu_index)] 
		pvalue=np.array(clu_p)[list(min_clu_index)] 
		node_pvalue.append(pvalue)
		min_path=clusters[index][min_clu_index[0]]
		pathway.append(path)
		color=np.repeat(colors[i],len(path))
		cluster=np.repeat((i+1),len(path)) 
		#print(i)
		#print(len(path))
		min_pathway.append(min_path)
		i=i+1
		node_color.append(color)
		cluster_id.append(cluster)

	min_pathway_name=list(np.array(name)[min_pathway]) 
	if database=="All":
		pass
	else:
		min_pathway_name = [c.replace('_PID','') for c in min_pathway_name]
		min_pathway_name = [c.replace('_KEGG','') for c in min_pathway_name]
		min_pathway_name = [c.replace('_Reactome','') for c in min_pathway_name]
		min_pathway_name = [c.replace('_INOH','') for c in min_pathway_name]
		min_pathway_name = [c.replace('_PANTHER','') for c in min_pathway_name]
		min_pathway_name = [c.replace('_NetPath','') for c in min_pathway_name]
		min_pathway_name = [c.replace('_HumanCyc','') for c in min_pathway_name]
		min_pathway_name = [c.replace('_WikiPathways','') for c in min_pathway_name]

	pathway_all = list(itertools.chain.from_iterable(pathway)) 
	cluster_id = list(itertools.chain.from_iterable(cluster_id))
	cluster_out=used.iloc[pathway_all,:]
	cluster_out['cluster']=cluster_id
	#cluster_out.to_csv(infofile,index=False,sep="\t")
	node_color = list(itertools.chain.from_iterable(node_color)) 
	node_pvalue = list(itertools.chain.from_iterable(node_pvalue)) 
	temp=new_matrix[pathway_all] 
	plot_matrix=temp[:,pathway_all] 
	if(np.all(plot_matrix == 0)):
		#plot_matrix=np.identity(len(pathway_all))
		file = open(outfile+".null",'w')
		file.close()
	else:
		f=open(topfile,'w') 
		for line in min_pathway_name:
			f.write(line+'\n')

		f.close()
		cluster_out.to_csv(infofile,index=False,sep="\t")
		plot_matrix=pd.DataFrame(plot_matrix)
		G = nx.from_pandas_adjacency(plot_matrix) 
		node_color_dict =dict(zip(range(0,len(node_color)),node_color))
		node_pvalue_dict =dict(zip(range(0,len(node_pvalue)),node_pvalue))
		nx.set_node_attributes(G, node_color_dict, 'nodecolor') 
		nx.set_node_attributes(G, node_pvalue_dict, 'pvalue')
		p4c.cytoscape_ping()
		p4c.get_network_list()
		p4c.delete_all_networks()
		p4c.command_sleep(2)
		n = p4c.create_network_from_networkx(G,title="mynetwork")
		time.sleep(4)
		try:
			p4c.delete_visual_style('CTpathway')
		except Exception as e:
			style_name = "CTpathway"
		else:
			style_name = "CTpathway"
		time.sleep(2)
		p4c.hide_all_panels()
		p4c.fit_content()
		defaults = {'NODE_SHAPE': "ELLIPSE", 'NODE_SIZE': 40, 'EDGE_TRANSPARENCY': 200,'EDGE_WIDTH': 5,'NODE_BORDER_PAINT':"#393B79",'NODE_BORDER_WIDTH':7,'NETWORK_SCALE_FACTOR':1.5}
		time.sleep(2)
		nodeFills = p4c.map_visual_property('node fill color','nodecolor','d',colors, colors)
		time.sleep(2)
		p4c.create_visual_style(style_name, defaults, [nodeFills])
		time.sleep(2)
		p4c.set_edge_color_default('#969696',style_name="CTpathway")
		time.sleep(2)
		p4c.set_visual_style(style_name) 
		time.sleep(2)
		p4c.bundle_edges(n)
		time.sleep(2)
		p4c.export_image(outfile, type='SVG', overwrite_file=True,network='mynetwork')
		p4c.delete_all_networks()
		p4c.cytoscape_free_memory()
		with open(outtemp, 'r') as f:
			all=f.read()

		f.close()
		buf = StringIO()
		buf.write(all)
		buf.seek(0,0) 
		out = StringIO()
		out.write(buf.readline())
		out.write(buf.readline())
		out.write(buf.readline())
		out.write(buf.readline())
		out.write(buf.readline())
		out.write(buf.readline())
		out.write(buf.readline())
		out.write(buf.readline())
		out.write(buf.readline())
		out.write(buf.readline())
		buf.readline() 
		buf.readline()
		b=''.join(buf.readlines()) 
		out.write(b)
		fd = open (outnetwork, 'w')
		fd.write(out.getvalue())
		fd.close()
		buf.close()
		out.close()
		fill_colors = ['#E41A1C', '#377EB8', '#4DAF4A', '#984EA3', '#FF7F00', '#FFFF33', '#A65628', '#F781BF', '#BD9E39', '#8DD3C7', '#FFFFB3', '#BEBADA', '#FB8072', '#80B1D3', '#FDB462', '#B3DE69', '#FCCDE5', '#FDD0A2', '#BC80BD', '#CCEBC5']
		data = pd.read_csv(topfile, sep='\t', header = None)
		pathway_len = []
		for i in data.iloc[:,0]:
			pathway_len.append(len(i))

		line = list()
		x = 10
		y = 10
		height = 22
		paper_width = 35 + 116*np.ceil(max(pathway_len)/10)
		paper_height = data.shape[0]*height + (data.shape[0]+1)*y
		paper_color = "white"
		line += svg_paper(paper_width,paper_height,color="white")
		for i in range(0, data.shape[0]):
			line += svg_rect("straight",x,y,height,height,fill_color=fill_colors[i])
			line += svg_txt(x+25,y+20,"20","black",data.iloc[i,0],textanchor=0)
			y += 30

		line += svg_end()
		f = open(outpathwayfile, "w+")
		for i in line:
			f.write(i)

		f.close()
		os.remove(topfile)
		os.remove(outtemp)
		#if index_use==0.3:
		#	out_pathway_pdf=dir+"/output/Clusters_in_Enrichment_map.pdf"
		#	cairosvg.svg2pdf(url=outpathwayfile,write_to=out_pathway_pdf)
		#	out_network_pdf=dir+"/output/Enrichment_map.pdf"
		#	networksvg=outfile+".svg"
		#	cairosvg.svg2pdf(url=networksvg,write_to=out_network_pdf)


def svg_paper(width,height,color="white"):
    import getpass
    import time
    svg_drawer = getpass.getuser()
    line = list()
    line.append("<?xml version=\"1.0\" encoding=\"iso-8859-1\"?>\n")
    line.append("<!DOCTYPE svg PUBLIC \"-//W3C//DTD SVG 20001102//EN\" \"http://www.w3.org/TR/2000/CR-SVG-20001102/DTD/svg-20001102.dtd\">\n\n")
    line.append("<svg viewBox=\"0 0 "+str(width)+" "+str(height)+"\">\n")
    line.append("<Drawer>"+svg_drawer+"</Drawer>\n")
    line.append("<Date>"+time.asctime(time.localtime(time.time()))+"</Date>\n")
    if color!="" :
        line.append("<rect x=\"0\" y=\"0\" width=\""+str(width)+"\" height=\""+str(height)+"\" fill=\""+color+"\"/>\n")
    return(line)

def svg_rect(lineType,x,y,width,height,fill_color,rx=1,ry=1,opacity=1):
    if lineType == "linear":
        line = ["<rect x=\""+str(x)+"\" y=\""+str(y)+"\" width=\""+str(width)+"\" height=\""+str(height)+"\" rx=\""+str(rx)+"\" ry=\""+str(ry)+"\" style=\"fill:"+fill_color+";opacity:"+str(opacity)+";\"/>\n"]
    elif lineType == "curve":
        line = ["<rect x=\""+str(x)+"\" y=\""+str(y)+"\" width=\""+str(width)+"\" height=\""+str(height)+"\" rx=\""+str(rx)+"\" ry=\""+str(ry)+"\" style=\"fill:"+fill_color+";opacity:"+str(opacity)+";\"/>\n"]
    elif lineType == "straight":
        line = ["<rect x=\""+str(x)+"\" y=\""+str(y)+"\" width=\""+str(width)+"\" height=\""+str(height)+"\" style=\"fill:"+fill_color+";opacity:"+str(opacity)+";\"/>\n"]
    return(line)

def svg_txt(x,y,size,color,text,vertical=0,textanchor=2,fontStyle="Normal"):
    if vertical == 0:
        svg_matrix = "1 0 0 1"
    elif vertical == 1:
        svg_matrix = "0 1 -1 0"
    elif vertical == 2:
        svg_matrix = "-1 0 0 -1"
    elif vertical == 3:
        svg_matrix = "0 -1 1 0"
    if textanchor == 0:
        anchor = "start"
    elif textanchor == 1:
        anchor = "middle"
    elif textanchor == 2:
        anchor = "end"
    line = ["<text fill=\""+color+"\" text-anchor=\""+anchor+"\" transform=\"matrix("+svg_matrix+" "+str(x)+" "+str(y)+")\" font-family=\"DejaVu Sans\" font-size=\""+str(size)+"\" font-weight=\""+str("500")+"\" font-style=\""+fontStyle+"\">"+str(text)+"</text>\n"]
    return(line)

def svg_end():
    return ["</svg>\n"]

def adjust_pvalues(pvalues, method='fdr'):
    return robjects.r['p.adjust'](robjects.FloatVector(pvalues), method=method)

###### Calculating gene risk score
###### Calculating pathway risk score
###### Identifying significantly risk pathways
def calculate_pathway_score_symbol(seed_file, database, sparse_matrix_file, random_number,data_type,sig_method,sig_value,dir):
	path= dir+"/output"
	outfile = path + "/pathway_result.txt"
	riskoutfile = path + "/risk_score_result.txt"
	riskgeneout = path + "/risk_gene_100.txt"
	sigfile = path + "/significant_pathway_result.txt"
	nodesfile = dir +'/data/big_network_nodes_symbol_name.txt'
	nodes=np.loadtxt(nodesfile,dtype='str')
	websitelink=dir+"/data/"+ database + "_website.txt"
	databasefile=dir+"/data/"+ database + "_pathway_index.txt"
	try:
		infile = open(seed_file,"r+")
		df =infile.read().splitlines()
		seed1=pd.DataFrame(df)
		data = pd.concat([seed1[0].str.split("\t", expand=True)], axis=1)
		seed1=data
		if(data.iloc[1,0].isdigit()):
			email_content = 'Sorry, your CTpathway job is failed! The Gene ID or Data Type in the input file does not match the selected parameter. Please check your input file and selected parameter!'
			print(email_content)
			exit(0)

		k_value=2
		if data_type == "Both":
			if seed1.shape[1] != 3:
				print(email_content)
				exit(0)

			seed1.columns=["n","fc","pvalue"]
			seed1['fc']=pd.to_numeric(seed1['fc']) #print(seed1.dtypes)
			seed1['pvalue']=pd.to_numeric(seed1['pvalue'])
			seed1['fc']=2**seed1.iloc[:,1]
			up=seed1[seed1.iloc[:,1]>1]
			up['new_fc']=(1-up['fc']**(-1/k_value))*(1-up['pvalue'])
			down=seed1[seed1.iloc[:,1]<1]
			down['new_fc']=(1-down['fc']**(1/k_value))*(1-down['pvalue'])
			ab1=pd.concat([up,down],axis=0)
		elif data_type == "Log2_FC":
			if seed1.shape[1] != 2:
				email_content = 'Sorry, your CTpathway job is failed! The Gene ID or Data Type in the input file does not match the selected parameter. Please check your input file and selected parameter!'
				print(email_content)
				exit(0)

			seed1.columns=["n","fc"]
			seed1['fc']=pd.to_numeric(seed1['fc'])
			seed1['fc']=2**seed1.iloc[:,1]
			up=seed1[seed1.iloc[:,1]>1]
			up['new_fc']=(1-up['fc']**(-1/k_value))
			down=seed1[seed1.iloc[:,1]<1]
			down['new_fc']=(1-down['fc']**(1/k_value))
			ab1=pd.concat([up,down],axis=0)
		else: #pvalue
			if seed1.shape[1] != 2:
				email_content = 'Sorry, your CTpathway job is failed! Please check your input file! The Data Type in the input file does not match the selected parameter. Please check your input file and selected parameter!'
				print(email_content)
				exit(0)

			seed1.columns=["n","pvalue"]
			seed1['pvalue']=pd.to_numeric(seed1['pvalue'])
			seed1['new_fc']=1-seed1['pvalue']
			ab1=seed1

		nodes1=pd.DataFrame(nodes,columns=list(('n')))
		weight=pd.merge(nodes1,ab1,how='left',on='n')
		weight=weight.fillna(0)
		random_fc=weight.new_fc
		random_number=int(random_number)
		p_threshold = 10/random_number
		for i in range(random_number):
			k=weight.new_fc.sample(15292)
			random_fc=np.concatenate((random_fc,k),axis=0)

		weight=weight.values
		sparse_matrix_csc =sparse.load_npz(sparse_matrix_file)
		matrix_csc=sparse_matrix_csc.toarray()
		random_fc_T=np.array(random_fc).reshape(random_number+1,15292).T
		all_score=np.dot(matrix_csc,random_fc_T)
		RS=all_score[:,0]
		RS_out_result=pd.DataFrame({"Gene symbol":nodes,"Risk score":RS}) 
		RS_out_result=RS_out_result.sort_values(['Risk score'],ascending=[False])
		RS_out_result.to_csv(riskoutfile,index=False,header=True,sep="\t")
		#risk_gene_100=RS_out_result.iloc[0:100,0]
		#a=risk_gene_100.values.reshape(20,5)
		#np.savetxt(riskgeneout,a,fmt="%s")
		#a=risk_gene_100.values.reshape(5,20)
		#np.savetxt(riskgeneout,a.T,fmt="%s")
		fExtremes = importr('fExtremes')
		args=np.zeros(2)
		args[1]=p_threshold
		Rcode='''cal_p<-function(data_rank,arr){
		numT=250;
		options(scipen=200)
		data_rank <- sort(data_rank,decreasing=TRUE);
		ob=arr[1]
		p_thre=arr[2]
		ori_p <- (1+length(which(data_rank>=ob)))/length(data_rank);
		#return(ori_p)
		if(ori_p>1){ori_p <- 1};
		if(ori_p>p_thre){return(ori_p)};
		t_score <- (data_rank[numT]+data_rank[numT+1])/2;
		F_gpd <- function(all_s){
		all_s <- sort(all_s,decreasing=TRUE);
		while(test_gpd(all_s)<0.05 & length(all_s)>100){
		all_s <- all_s[1:(length(all_s)-10)];
		}
		fit_d <-  gpdFit(all_s,u=0);
		return(fit_d);
		}
		test_gpd <- function(all_s){
		fit_d <-  gpdFit(all_s,u=0);
		sim_d <-  as.numeric(gpdSim(model = list(beta=fit_d@fit$par.ests["beta"],mu=0,xi=fit_d@fit$par.ests["xi"]), 
		n = length(all_s)));
		pv<-ks.test(all_s,sim_d)$p.value;return(pv);
		} 
		fit_d <- try(F_gpd(c(ob,data_rank[1:numT])-t_score));
		if(class(fit_d)=="try-error"){return(ori_p);}
		pv <- pgpd(ob-t_score,mu=0,beta=fit_d@fit$par.ests["beta"],xi=fit_d@fit$par.ests["xi"]);
		adj_p <- ori_p*(1-pv[1]);
		return(as.numeric(adj_p));
		}'''
		pareto=r(Rcode)
		pathway_name=list()
		p_value_vector=list()
		score_vector=list()
		node_number=list()
		for line in open(databasefile):
			line=line.strip('\n').split('\t')
			pathway_name.append(line[0])
			score_total=np.mean(all_score[list(map(int, line[1:]))],0)
			score_vector.append(score_total[0])
			p_value=len(score_total[score_total>score_total[0]])/random_number
			if p_value<=p_threshold:
				args[0]=score_total[0]
				p_value=np.array(pareto(robjects.FloatVector(score_total[1:(random_number+1)]),robjects.FloatVector(args)))[0]

			
			p_value_vector.append(p_value)
			node_number.append(len(line)-1)

		out_result=pd.DataFrame({"pathway_name":pathway_name,"PS":score_vector,"p_value":p_value_vector,"node_number":node_number})
		out_result=out_result.sort_values(['p_value','PS'],ascending=[True,False])
		adj_p=adjust_pvalues(out_result.p_value)
		out_result['FDR']=adj_p
		out_result['Source']=database
		order=['pathway_name','Source','PS','p_value','FDR','node_number']
		out_result=out_result[order]
		out_result['p_value']=out_result['p_value'].map(lambda x:('%.3e')%x)
		out_result['FDR']=out_result['FDR'].map(lambda x:('%.3e')%x)
		out_result['PS']=out_result['PS'].map(lambda x:('%.4f')%x)
		link=pd.read_csv(websitelink,sep="\t")
		out_result=pd.merge(out_result,link,how='left',on='pathway_name')
		out_result=out_result.fillna("null")
		out_result.to_csv(outfile,index=False,header=True,sep="\t")
		if sig_method=="FDR":
			used= out_result.loc[out_result['FDR'].astype('float')< float(sig_value)]
		else:
			used= out_result.loc[out_result['p_value'].astype('float')< float(sig_value)]

		used.to_csv(sigfile,index=False,header=True,sep="\t")
		path_svg= path
		out_svg_1 = path_svg + "/bar_20.svg"
		out_png_1 = path + "/bar.pdf"
		out_bubble_svg_1 = path_svg + "/bubble_20.svg"
		out_bubble_png_1 = path + "/bubble.pdf"
		plt_bubble_bar_svg(outfile,out_bubble_svg_1,out_svg_1,20,"single")
		cairosvg.svg2pdf(url=out_bubble_svg_1,write_to=out_bubble_png_1)
		cairosvg.svg2pdf(url=out_svg_1,write_to=out_png_1)
		os.remove(out_bubble_svg_1)
		os.remove(out_svg_1)
		index=0.3
		index_list=[0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9]
		for i in index_list:
			cytoscape(i,sig_method,sig_value,"single",dir)
			cytoscape_file=Path(path_svg + "/network_temp_"+str(i)+".null")
			if cytoscape_file.is_file():
				os.remove(cytoscape_file)
				break


		email_content ="Congratulations! Your CTpathway job has been completed."
		print(email_content)
	except Exception as e:
		error_log=path + '/error.txt'
		traceback.print_exc(file=open(error_log,'a+'))
		email_content = 'Sorry, your CTpathway job is failed! Please check your input file!'
		print(email_content)
		#shutil.rmtree(path)

def calculate_pathway_score_symbol_all(seed_file, database, sparse_matrix_file, random_number,data_type,sig_method,sig_value,dir):
	path= dir+"/output"
	outfile = path + "/pathway_result.txt"
	riskoutfile = path + "/risk_score_result.txt"
	riskgeneout = path + "/risk_gene_100.txt"
	sigfile = path + "/significant_pathway_result.txt"
	nodesfile = dir +'/data/big_network_nodes_symbol_name.txt'
	nodes=np.loadtxt(nodesfile,dtype='str')
	try:
		infile = open(seed_file,"r+")
		#df = infile.readlines()
		df =infile.read().splitlines()
		seed1=pd.DataFrame(df)
		data = pd.concat([seed1[0].str.split("\t", expand=True)], axis=1)
		seed1=data
		if(data.iloc[1,0].isdigit()):
			email_content = 'Sorry, your CTpathway job is failed! The Gene ID or Data Type in the input file does not match the selected parameter. Please check your input file and selected parameter!'
			print(email_content)
			exit(0)

		k_value=2
		if data_type == "Both":
			if seed1.shape[1] != 3:
				email_content = 'Sorry, your CTpathway job is failed! The Gene ID or Data Type in the input file does not match the selected parameter. Please check your input file and selected parameter!'
				print(email_content)
				exit(0)

			seed1.columns=["n","fc","pvalue"]
			seed1['fc']=pd.to_numeric(seed1['fc']) #print(seed1.dtypes)
			seed1['pvalue']=pd.to_numeric(seed1['pvalue'])
			seed1['fc']=2**seed1.iloc[:,1]
			up=seed1[seed1.iloc[:,1]>1]
			up['new_fc']=(1-up['fc']**(-1/k_value))*(1-up['pvalue'])
			down=seed1[seed1.iloc[:,1]<1]
			down['new_fc']=(1-down['fc']**(1/k_value))*(1-down['pvalue'])
			ab1=pd.concat([up,down],axis=0)
		elif data_type == "Log2_FC":
			if seed1.shape[1] != 2:
				email_content = 'Sorry, your CTpathway job is failed! The Gene ID or Data Type in the input file does not match the selected parameter. Please check your input file and selected parameter!'
				print(email_content)
				exit(0)

			seed1.columns=["n","fc"]
			seed1['fc']=pd.to_numeric(seed1['fc'])
			seed1['fc']=2**seed1.iloc[:,1]
			up=seed1[seed1.iloc[:,1]>1]
			up['new_fc']=(1-up['fc']**(-1/k_value))
			down=seed1[seed1.iloc[:,1]<1]
			down['new_fc']=(1-down['fc']**(1/k_value))
			ab1=pd.concat([up,down],axis=0)
		else: #pvalue
			if seed1.shape[1] != 2:
				email_content = 'Sorry, your CTpathway job is failed! The Gene ID or Data Type in the input file does not match the selected parameter. Please check your input file and selected parameter!'
				print(email_content)
				exit(0)

			seed1.columns=["n","pvalue"]
			seed1['pvalue']=pd.to_numeric(seed1['pvalue'])
			seed1['new_fc']=1-seed1['pvalue']
			ab1=seed1

		nodes1=pd.DataFrame(nodes,columns=list(('n')))
		weight=pd.merge(nodes1,ab1,how='left',on='n')
		weight=weight.fillna(0)
		random_fc=weight.new_fc
		random_number=int(random_number)
		p_threshold = 10/random_number
		for i in range(random_number):
			k=weight.new_fc.sample(15292)
			random_fc=np.concatenate((random_fc,k),axis=0)

		weight=weight.values
		sparse_matrix_csc =sparse.load_npz(sparse_matrix_file)
		matrix_csc=sparse_matrix_csc.toarray()
		random_fc_T=np.array(random_fc).reshape(random_number+1,15292).T
		all_score=np.dot(matrix_csc,random_fc_T)
		RS=all_score[:,0]
		RS_out_result=pd.DataFrame({"Gene symbol":nodes,"Risk score":RS}) 
		RS_out_result=RS_out_result.sort_values(['Risk score'],ascending=[False])
		RS_out_result.to_csv(riskoutfile,index=False,header=True,sep="\t")
		#risk_gene_100=RS_out_result.iloc[0:100,0]
		#a=risk_gene_100.values.reshape(20,5)
		#np.savetxt(riskgeneout,a,fmt="%s")
		#a=risk_gene_100.values.reshape(5,20)
		#np.savetxt(riskgeneout,a.T,fmt="%s")
		fExtremes = importr('fExtremes')
		args=np.zeros(2)
		args[1]=p_threshold
		Rcode='''cal_p<-function(data_rank,arr){
		numT=250;
		options(scipen=200)
		data_rank <- sort(data_rank,decreasing=TRUE);
		ob=arr[1]
		p_thre=arr[2]
		ori_p <- (1+length(which(data_rank>=ob)))/length(data_rank);
		#return(ori_p)
		if(ori_p>1){ori_p <- 1};
		if(ori_p>p_thre){return(ori_p)};
		t_score <- (data_rank[numT]+data_rank[numT+1])/2;
		F_gpd <- function(all_s){
		all_s <- sort(all_s,decreasing=TRUE);
		while(test_gpd(all_s)<0.05 & length(all_s)>100){
		all_s <- all_s[1:(length(all_s)-10)];
		}
		fit_d <-  gpdFit(all_s,u=0);
		return(fit_d);
		}
		test_gpd <- function(all_s){
		fit_d <-  gpdFit(all_s,u=0);
		sim_d <-  as.numeric(gpdSim(model = list(beta=fit_d@fit$par.ests["beta"],mu=0,xi=fit_d@fit$par.ests["xi"]), 
		n = length(all_s)));
		pv<-ks.test(all_s,sim_d)$p.value;return(pv);
		} 
		fit_d <- try(F_gpd(c(ob,data_rank[1:numT])-t_score));
		if(class(fit_d)=="try-error"){return(ori_p);}
		pv <- pgpd(ob-t_score,mu=0,beta=fit_d@fit$par.ests["beta"],xi=fit_d@fit$par.ests["xi"]);
		adj_p <- ori_p*(1-pv[1]);
		return(as.numeric(adj_p));
		}'''
		pareto=r(Rcode)
###############KEGG
		pathway_name=list()
		p_value_vector=list()
		score_vector=list()
		node_number=list()
		databasefile=dir+"/data/"+ "KEGG_pathway_index.txt"
		for line in open(databasefile):
			line=line.strip('\n').split('\t')
			pathway_name.append(line[0])
			score_total=np.mean(all_score[list(map(int, line[1:]))],0)
			score_vector.append(score_total[0])
			p_value=len(score_total[score_total>score_total[0]])/random_number
			if p_value<=p_threshold:
				args[0]=score_total[0]
				p_value=np.array(pareto(robjects.FloatVector(score_total[1:(random_number+1)]),robjects.FloatVector(args)))[0]

			p_value_vector.append(p_value)
			node_number.append(len(line)-1)

		kegg=pd.DataFrame({"pathway_name":pathway_name,"PS":score_vector,"p_value":p_value_vector,"node_number":node_number})
		kegg['Source']="KEGG"
		websitelink=dir+"/data/KEGG_website.txt"
		link=pd.read_csv(websitelink,sep="\t")
		kegg=pd.merge(kegg,link,how='left',on='pathway_name')
		kegg=kegg.fillna("null")
###############Reactome
		pathway_name=list()
		p_value_vector=list()
		score_vector=list()
		node_number=list()
		databasefile=dir+"/data/"+ "Reactome_pathway_index.txt"
		for line in open(databasefile):
			line=line.strip('\n').split('\t')
			pathway_name.append(line[0])
			score_total=np.mean(all_score[list(map(int, line[1:]))],0)
			score_vector.append(score_total[0])
			p_value=len(score_total[score_total>score_total[0]])/random_number
			if p_value<=p_threshold:
				args[0]=score_total[0]
				p_value=np.array(pareto(robjects.FloatVector(score_total[1:(random_number+1)]),robjects.FloatVector(args)))[0]

			p_value_vector.append(p_value)
			node_number.append(len(line)-1)

		Reactome=pd.DataFrame({"pathway_name":pathway_name,"PS":score_vector,"p_value":p_value_vector,"node_number":node_number})
		Reactome['Source']="Reactome"
		websitelink=dir+"/data/Reactome_website.txt"
		link=pd.read_csv(websitelink,sep="\t")
		Reactome=pd.merge(Reactome,link,how='left',on='pathway_name')
		Reactome=Reactome.fillna("null")
###########Pid
		pathway_name=list()
		p_value_vector=list()
		score_vector=list()
		node_number=list()
		databasefile=dir+"/data/"+ "PID_pathway_index.txt"
		for line in open(databasefile):
			line=line.strip('\n').split('\t')
			pathway_name.append(line[0])
			score_total=np.mean(all_score[list(map(int, line[1:]))],0)
			score_vector.append(score_total[0])
			p_value=len(score_total[score_total>score_total[0]])/random_number
			if p_value<=p_threshold:
				args[0]=score_total[0]
				p_value=np.array(pareto(robjects.FloatVector(score_total[1:(random_number+1)]),robjects.FloatVector(args)))[0]
            
			p_value_vector.append(p_value)
			node_number.append(len(line)-1)

		Pid=pd.DataFrame({"pathway_name":pathway_name,"PS":score_vector,"p_value":p_value_vector,"node_number":node_number})
		Pid['Source']="PID"
		websitelink=dir+"/data/PID_website.txt"
		link=pd.read_csv(websitelink,sep="\t")
		Pid=pd.merge(Pid,link,how='left',on='pathway_name')
		Pid=Pid.fillna("null")
###########Netpath
		pathway_name=list()
		p_value_vector=list()
		score_vector=list()
		node_number=list()
		databasefile=dir+"/data/"+ "NetPath_pathway_index.txt"
		for line in open(databasefile):
			line=line.strip('\n').split('\t')
			pathway_name.append(line[0])
			score_total=np.mean(all_score[list(map(int, line[1:]))],0)
			score_vector.append(score_total[0])
			p_value=len(score_total[score_total>score_total[0]])/random_number
			if p_value<=p_threshold:
				args[0]=score_total[0]
				p_value=np.array(pareto(robjects.FloatVector(score_total[1:(random_number+1)]),robjects.FloatVector(args)))[0]

			p_value_vector.append(p_value)
			node_number.append(len(line)-1)

		Netpath=pd.DataFrame({"pathway_name":pathway_name,"PS":score_vector,"p_value":p_value_vector,"node_number":node_number})
		Netpath['Source']="NetPath"
		websitelink=dir+"/data/NetPath_website.txt"
		link=pd.read_csv(websitelink,sep="\t")
		Netpath=pd.merge(Netpath,link,how='left',on='pathway_name')
		Netpath=Netpath.fillna("null")
#########Inoh
		pathway_name=list()
		p_value_vector=list()
		score_vector=list()
		node_number=list()
		databasefile=dir+"/data/"+ "INOH_pathway_index.txt"
		for line in open(databasefile):
			line=line.strip('\n').split('\t')
			pathway_name.append(line[0])
			score_total=np.mean(all_score[list(map(int, line[1:]))],0)
			score_vector.append(score_total[0])
			p_value=len(score_total[score_total>score_total[0]])/random_number
			if p_value<=p_threshold:
				args[0]=score_total[0]
				p_value=np.array(pareto(robjects.FloatVector(score_total[1:(random_number+1)]),robjects.FloatVector(args)))[0]

			p_value_vector.append(p_value)
			node_number.append(len(line)-1)

		Inoh=pd.DataFrame({"pathway_name":pathway_name,"PS":score_vector,"p_value":p_value_vector,"node_number":node_number})
		Inoh['Source']="INOH"
		websitelink=dir+"/data/INOH_website.txt"
		link=pd.read_csv(websitelink,sep="\t")
		Inoh=pd.merge(Inoh,link,how='left',on='pathway_name')
		Inoh=Inoh.fillna("null")
#########Humancyc
		pathway_name=list()
		p_value_vector=list()
		score_vector=list()
		node_number=list()
		databasefile=dir+"/data/"+ "HumanCyc_pathway_index.txt"
		for line in open(databasefile):
			line=line.strip('\n').split('\t')
			pathway_name.append(line[0])
			score_total=np.mean(all_score[list(map(int, line[1:]))],0)
			score_vector.append(score_total[0])
			p_value=len(score_total[score_total>score_total[0]])/random_number
			if p_value<=p_threshold:
				args[0]=score_total[0]
				p_value=np.array(pareto(robjects.FloatVector(score_total[1:(random_number+1)]),robjects.FloatVector(args)))[0]

			p_value_vector.append(p_value)
			node_number.append(len(line)-1)

		Humancyc=pd.DataFrame({"pathway_name":pathway_name,"PS":score_vector,"p_value":p_value_vector,"node_number":node_number})
		Humancyc['Source']="HumanCyc"
		websitelink=dir+"/data/HumanCyc_website.txt"
		link=pd.read_csv(websitelink,sep="\t")
		Humancyc=pd.merge(Humancyc,link,how='left',on='pathway_name')
		Humancyc=Humancyc.fillna("null")
##########Panther
		pathway_name=list()
		p_value_vector=list()
		score_vector=list()
		node_number=list()
		databasefile=dir+"/data/"+ "PANTHER_pathway_index.txt"
		for line in open(databasefile):
			line=line.strip('\n').split('\t')
			pathway_name.append(line[0])
			score_total=np.mean(all_score[list(map(int, line[1:]))],0)
			score_vector.append(score_total[0])
			p_value=len(score_total[score_total>score_total[0]])/random_number
			if p_value<=p_threshold:
				args[0]=score_total[0]
				p_value=np.array(pareto(robjects.FloatVector(score_total[1:(random_number+1)]),robjects.FloatVector(args)))[0]

			p_value_vector.append(p_value)
			node_number.append(len(line)-1)

		Panther=pd.DataFrame({"pathway_name":pathway_name,"PS":score_vector,"p_value":p_value_vector,"node_number":node_number})
		Panther['Source']="PANTHER"
		websitelink=dir+"/data/PANTHER_website.txt"
		link=pd.read_csv(websitelink,sep="\t")
		Panther=pd.merge(Panther,link,how='left',on='pathway_name')
		Panther=Panther.fillna("null")
##########Wikipathways
		pathway_name=list()
		p_value_vector=list()
		score_vector=list()
		node_number=list()
		databasefile=dir+"/data/"+ "WikiPathways_pathway_index.txt"
		for line in open(databasefile):
			line=line.strip('\n').split('\t')
			pathway_name.append(line[0])
			score_total=np.mean(all_score[list(map(int, line[1:]))],0)
			score_vector.append(score_total[0])
			p_value=len(score_total[score_total>score_total[0]])/random_number
			if p_value<=p_threshold:
				args[0]=score_total[0]
				p_value=np.array(pareto(robjects.FloatVector(score_total[1:(random_number+1)]),robjects.FloatVector(args)))[0]

			p_value_vector.append(p_value)
			node_number.append(len(line)-1)

		Wikipathways=pd.DataFrame({"pathway_name":pathway_name,"PS":score_vector,"p_value":p_value_vector,"node_number":node_number})
		Wikipathways['Source']="WikiPathways"
		websitelink=dir+"/data/WikiPathways_website.txt"
		link=pd.read_csv(websitelink,sep="\t")
		Wikipathways=pd.merge(Wikipathways,link,how='left',on='pathway_name')
		Wikipathways=Wikipathways.fillna("null")
		out_result=pd.concat([kegg,Inoh,Reactome,Netpath,Pid,Panther,Wikipathways,Humancyc])
		out_result=out_result.sort_values(['p_value','PS'],ascending=[True,False])
		adj_p=adjust_pvalues(out_result.p_value)
		out_result['FDR']=adj_p
		order=['pathway_name','Source','PS','p_value','FDR','node_number',"website"]
		out_result=out_result[order]
		out_result['p_value']=out_result['p_value'].map(lambda x:('%.3e')%x)
		out_result['FDR']=out_result['FDR'].map(lambda x:('%.3e')%x)
		out_result['PS']=out_result['PS'].map(lambda x:('%.4f')%x)
		out_result=out_result[order]
		out_result.to_csv(outfile,index=False,header=True,sep="\t")
		if sig_method=="FDR":
			used= out_result.loc[out_result['FDR'].astype('float')< float(sig_value)]
		else:
			used= out_result.loc[out_result['p_value'].astype('float')< float(sig_value)]

		used.to_csv(sigfile,index=False,header=True,sep="\t")
		path_svg= path
		out_svg_1 = path_svg + "/bar_20.svg"
		out_png_1 = path + "/bar.pdf"
		out_bubble_svg_1 = path_svg + "/bubble_20.svg"
		out_bubble_png_1 = path + "/bubble.pdf"
		plt_bubble_bar_svg(outfile,out_bubble_svg_1,out_svg_1,20,"All")
		cairosvg.svg2pdf(url=out_bubble_svg_1,write_to=out_bubble_png_1)
		cairosvg.svg2pdf(url=out_svg_1,write_to=out_png_1)
		os.remove(out_bubble_svg_1)
		os.remove(out_svg_1)
		index=0.3
		index_list=[0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9]
		for i in index_list:
			cytoscape(i,sig_method,sig_value,"All",dir)
			cytoscape_file=Path(path_svg + "/network_temp_"+str(i)+".null")
			if cytoscape_file.is_file():
				os.remove(cytoscape_file)
				break


		email_content ='Congratulations! Your CTpathway job has been completed. '
		print(email_content)
	except Exception as e:
		error_log=path + '/error.txt'
		traceback.print_exc(file=open(error_log,'a+'))
		email_content = 'Sorry, your CTpathway job is failed! Please check your input file!'
		print(email_content)
		#shutil.rmtree(path)

def calculate_pathway_score(seed_file, database, sparse_matrix_file, random_number,data_type,sig_method,sig_value,dir):
	path= dir+"/output"
	outfile = path + "/pathway_result.txt"
	riskoutfile = path + "/risk_score_result.txt"
	riskgeneout = path + "/risk_gene_100.txt"
	sigfile = path + "/significant_pathway_result.txt"
	nodesfile = dir +'/data/big_network_nodes.txt'
	nodes=np.loadtxt(nodesfile)
	websitelink=dir+"/data/"+ database + "_website.txt"
	databasefile=dir+"/data/"+ database + "_pathway_index.txt"
	try:
		seed=np.loadtxt(seed_file)
	except Exception as e:
		error_log=path + '/error.txt'
		traceback.print_exc(file=open(error_log,'a+'))
		email_content = 'Sorry, your CTpathway job is failed! The Gene ID or Data Type in the input file does not match the selected parameter. Please check your input file and selected parameter!'
		print(email_content)
		exit()

	try:
		#seed=np.loadtxt(seed_file) 
		k_value = 2
		if data_type == "Both":
			if seed.shape[1] != 3:
				print(email_content)
				exit(0)

			seed[:,1]=2**seed[:,1]
			seed=seed[:,[0,1,2]]
			up=seed[seed[:,1]>=1,:]
			up[:,1]=(1-up[:,1]**(-1/k_value))*(1-up[:,2])
			down=seed[seed[:,1]<1,:]
			down[:,1]=(1-down[:,1]**(1/k_value))*(1-down[:,2])
			ab=np.concatenate((up, down), axis=0)
			ab1 = pd.DataFrame(ab, columns=list(('n','fc','pvalue')))
		elif data_type == "Log2_FC":
			if seed.shape[1] != 2:
				email_content = 'Sorry, your CTpathway job is failed! The Gene ID or Data Type in the input file does not match the selected parameter. Please check your input file and selected parameter!'
				print(email_content)
				exit(0)

			seed[:,1]=2**seed[:,1]
			seed=seed[:,[0,1]]
			up=seed[seed[:,1]>=1,:]
			up[:,1]=(1-up[:,1]**(-1/k_value))
			down=seed[seed[:,1]<1,:]
			down[:,1]=(1-down[:,1]**(1/k_value))
			ab=np.concatenate((up, down), axis=0)
			ab1 = pd.DataFrame(ab, columns=list(('n','fc')))
		else:
			if seed.shape[1] != 2:
				email_content = 'Sorry, your CTpathway job is failed! Please check your input file! The Data Type in the input file does not match the selected parameter. Please check your input file and selected parameter!'
				print(email_content)
				exit(0)

			seed[:,1]=1-seed[:,1]
			ab1 = pd.DataFrame(seed, columns=list(('n','fc')))

		nodes1 = pd.DataFrame(nodes, columns=list(('n')))
		weight=pd.merge(nodes1, ab1, how='left', on='n')
		weight=weight.fillna(0)
		random_fc=weight.fc
		random_number=int(random_number)
		p_threshold = 10/random_number
		for i in range(random_number):
			k=weight.fc.sample(15292)
			random_fc=np.concatenate((random_fc, k), axis=0)

		sparse_matrix_csc =sparse.load_npz(sparse_matrix_file)
		matrix_csc=sparse_matrix_csc.toarray()
		random_fc_T=np.array(random_fc).reshape(random_number+1, 15292).T
		all_score=np.dot(matrix_csc, random_fc_T)
		RS=all_score[:,0]
		RS_out_result=pd.DataFrame({"Entrez ID":nodes.astype(np.int),"Risk score":RS})
		RS_out_result=RS_out_result.sort_values(['Risk score'],ascending=[False])
		RS_out_result.to_csv(riskoutfile,index=False,header=True,sep="\t")
		#risk_gene_100=RS_out_result.iloc[0:100,0]
		#a=risk_gene_100.values.reshape(20,5)
		#np.savetxt(riskgeneout,a,fmt="%s")
		#a=risk_gene_100.values.reshape(5,20)
		#np.savetxt(riskgeneout,a.T,fmt="%s")
		fExtremes = importr('fExtremes')
		args=np.zeros(2)
		args[1]=p_threshold
		Rcode='''cal_p<-function(data_rank,arr){
		numT=250;
		options(scipen=200)
		data_rank <- sort(data_rank,decreasing=TRUE);
		ob=arr[1]
		p_thre=arr[2]
		ori_p <- (1+length(which(data_rank>=ob)))/length(data_rank);
		#return(ori_p)
		if(ori_p>1){ori_p <- 1};
		if(ori_p>p_thre){return(ori_p)};
		t_score <- (data_rank[numT]+data_rank[numT+1])/2;
		F_gpd <- function(all_s){
		all_s <- sort(all_s,decreasing=TRUE);
		while(test_gpd(all_s)<0.05 & length(all_s)>100){
		all_s <- all_s[1:(length(all_s)-10)];
		}
		fit_d <-  gpdFit(all_s,u=0);
		return(fit_d);
		}
		test_gpd <- function(all_s){
		fit_d <-  gpdFit(all_s,u=0);
		sim_d <-  as.numeric(gpdSim(model = list(beta=fit_d@fit$par.ests["beta"],mu=0,xi=fit_d@fit$par.ests["xi"]),
		n = length(all_s)));
		pv<-ks.test(all_s,sim_d)$p.value;return(pv);
		}
		fit_d <- try(F_gpd(c(ob,data_rank[1:numT])-t_score));
		if(class(fit_d)=="try-error"){return(ori_p);}
		pv <- pgpd(ob-t_score,mu=0,beta=fit_d@fit$par.ests["beta"],xi=fit_d@fit$par.ests["xi"]);
		adj_p <- ori_p*(1-pv[1]);
		return(as.numeric(adj_p));
		}'''
		pareto=r(Rcode)
		pathway_name=list()
		p_value_vector=list()
		score_vector=list()
		node_number=list()
		for line in open(databasefile):
			line=line.strip('\n').split('\t')
			pathway_name.append(line[0])
			score_total=np.mean(all_score[list(map(int, line[1:]))],0)
			score_vector.append(score_total[0])
			p_value=len(score_total[score_total>score_total[0]])/random_number
			if p_value<=p_threshold:
				args[0]=score_total[0]
				p_value=np.array(pareto(robjects.FloatVector(score_total[1:(random_number+1)]),robjects.FloatVector(args)))[0]

			p_value_vector.append(p_value)
			node_number.append(len(line)-1)

		out_result=pd.DataFrame({"pathway_name":pathway_name,"PS":score_vector,"p_value":p_value_vector,"node_number":node_number})
		out_result=out_result.sort_values(['p_value','PS'],ascending=[True,False])
		adj_p=adjust_pvalues(out_result.p_value)
		out_result['FDR']=adj_p
		out_result['Source']=database
		order=['pathway_name','Source','PS','p_value','FDR','node_number']
		out_result=out_result[order]
		out_result['p_value']=out_result['p_value'].map(lambda x:('%.3e')%x)
		out_result['FDR']=out_result['FDR'].map(lambda x:('%.3e')%x)
		out_result['PS']=out_result['PS'].map(lambda x:('%.4f')%x)
		out_result=out_result[order]
		link=pd.read_csv(websitelink,sep="\t")
		out_result=pd.merge(out_result,link,how='left',on='pathway_name')
		out_result=out_result.fillna("null")
		out_result.to_csv(outfile,index=False,header=True,sep="\t")
		if sig_method=="FDR":
			used= out_result.loc[out_result['FDR'].astype('float')< float(sig_value)]
		else:
			used= out_result.loc[out_result['p_value'].astype('float')< float(sig_value)]

		used.to_csv(sigfile,index=False,header=True,sep="\t")
		path_svg= path
		out_svg_1 = path_svg + "/bar_20.svg"
		out_png_1 = path + "/bar.pdf"
		out_bubble_svg_1 = path_svg + "/bubble_20.svg"
		out_bubble_png_1 = path + "/bubble.pdf"
		plt_bubble_bar_svg(outfile,out_bubble_svg_1,out_svg_1,20,"single")
		cairosvg.svg2pdf(url=out_bubble_svg_1,write_to=out_bubble_png_1)
		cairosvg.svg2pdf(url=out_svg_1,write_to=out_png_1)
		os.remove(out_bubble_svg_1)
		os.remove(out_svg_1)
		index=0.3
		index_list=[0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9]
		for i in index_list:
			cytoscape(i,sig_method,sig_value,"single",dir)
			cytoscape_file=Path(path_svg + "/network_temp_"+str(i)+".null")
			if cytoscape_file.is_file():
				os.remove(cytoscape_file)
				break


		email_content ="Congratulations! Your CTpathway job has been completed."
		print(email_content)
	except Exception as e:
		error_log=path + '/error.txt'
		traceback.print_exc(file=open(error_log,'a+'))
		email_content = 'Sorry, your CTpathway job is failed! Please check your input file!'
		print(email_content)
		#shutil.rmtree(path)

def calculate_pathway_score_all(seed_file, database, sparse_matrix_file, random_number,data_type,sig_method,sig_value,dir):
	path= dir+"/output"
	outfile = path + "/pathway_result.txt"
	riskoutfile = path + "/risk_score_result.txt"
	riskgeneout = path + "/risk_gene_100.txt"
	sigfile = path + "/significant_pathway_result.txt"
	nodesfile = dir +'/data/big_network_nodes.txt'
	nodes=np.loadtxt(nodesfile)
	websitelink=dir+"/data/"+ database + "_website.txt"
	databasefile=dir+"/data/"+ database + "_pathway_index.txt"
	try:
		seed=np.loadtxt(seed_file)
	except Exception as e:
		error_log=path + '/error.txt'
		traceback.print_exc(file=open(error_log,'a+'))
		email_content = 'Sorry, your CTpathway job is failed! The Gene ID or Data Type in the input file does not match the selected parameter. Please check your input file and selected parameter!'
		print(email_content)
		exit()

	try:
		#seed=np.loadtxt(seed_file)
		k_value = 2
		if data_type == "Both":
			if seed.shape[1] != 3:
				email_content = 'Sorry, your CTpathway job is failed! The Gene ID or Data Type in the input file does not match the selected parameter. Please check your input file and selected parameter!'
				print(email_content)
				exit(0)

			seed[:,1]=2**seed[:,1]
			seed=seed[:,[0,1,2]]
			up=seed[seed[:,1]>=1,:]
			up[:,1]=(1-up[:,1]**(-1/k_value))*(1-up[:,2])
			down=seed[seed[:,1]<1,:]
			down[:,1]=(1-down[:,1]**(1/k_value))*(1-down[:,2])
			ab=np.concatenate((up, down), axis=0)
			ab1 = pd.DataFrame(ab, columns=list(('n','fc','pvalue')))
		elif data_type == "Log2_FC":
			if seed.shape[1] != 2:
				email_content = 'Sorry, your CTpathway job is failed! The Gene ID or Data Type in the input file does not match the selected parameter. Please check your input file and selected parameter!'
				print(email_content)
				exit(0)

			seed[:,1]=2**seed[:,1]
			seed=seed[:,[0,1]]
			up=seed[seed[:,1]>=1,:]
			up[:,1]=(1-up[:,1]**(-1/k_value))
			down=seed[seed[:,1]<1,:]
			down[:,1]=(1-down[:,1]**(1/k_value))
			ab=np.concatenate((up, down), axis=0)
			ab1 = pd.DataFrame(ab, columns=list(('n','fc')))
		else:
			if seed.shape[1] != 2:
				email_content = 'Sorry, your CTpathway job is failed! The Gene ID or Data Type in the input file does not match the selected parameter. Please check your input file and selected parameter!'
				print(email_content)
				exit(0)

			seed[:,1]=1-seed[:,1]
			ab1 = pd.DataFrame(seed, columns=list(('n','fc')))

		nodes1=pd.DataFrame(nodes, columns=list(('n')))
		weight=pd.merge(nodes1, ab1, how='left', on='n')
		weight=weight.fillna(0)
		random_fc=weight.fc
		random_number=int(random_number)
		p_threshold = 10/random_number
		for i in range(random_number):
			k=weight.fc.sample(15292)
			random_fc=np.concatenate((random_fc, k), axis=0)

		sparse_matrix_csc =sparse.load_npz(sparse_matrix_file)
		matrix_csc=sparse_matrix_csc.toarray()
		random_fc_T=np.array(random_fc).reshape(random_number+1, 15292).T
		all_score=np.dot(matrix_csc, random_fc_T)
		RS=all_score[:,0]
		RS_out_result=pd.DataFrame({"Entrez ID":nodes.astype(np.int),"Risk score":RS}) 
		RS_out_result=RS_out_result.sort_values(['Risk score'],ascending=[False])
		RS_out_result.to_csv(riskoutfile,index=False,header=True,sep="\t")
		#risk_gene_100=RS_out_result.iloc[0:100,0]
		#a=risk_gene_100.values.reshape(20,5)
		#np.savetxt(riskgeneout,a,fmt="%s")
		#a=risk_gene_100.values.reshape(5,20)
		#np.savetxt(riskgeneout,a.T,fmt="%s")
		fExtremes = importr('fExtremes')
		args=np.zeros(2)
		args[1]=p_threshold
		Rcode='''cal_p<-function(data_rank,arr){
		numT=250;
		options(scipen=200)
		data_rank <- sort(data_rank,decreasing=TRUE);
		ob=arr[1]
		p_thre=arr[2]
		ori_p <- (1+length(which(data_rank>=ob)))/length(data_rank);
		#return(ori_p)
		if(ori_p>1){ori_p <- 1};
		if(ori_p>p_thre){return(ori_p)};
		t_score <- (data_rank[numT]+data_rank[numT+1])/2;
		F_gpd <- function(all_s){
		all_s <- sort(all_s,decreasing=TRUE);
		while(test_gpd(all_s)<0.05 & length(all_s)>100){
		all_s <- all_s[1:(length(all_s)-10)];
		}
		fit_d <-  gpdFit(all_s,u=0);
		return(fit_d);
		}
		test_gpd <- function(all_s){
		fit_d <-  gpdFit(all_s,u=0);
		sim_d <-  as.numeric(gpdSim(model = list(beta=fit_d@fit$par.ests["beta"],mu=0,xi=fit_d@fit$par.ests["xi"]),
		n = length(all_s)));
		pv<-ks.test(all_s,sim_d)$p.value;return(pv);
		}
		fit_d <- try(F_gpd(c(ob,data_rank[1:numT])-t_score));
		if(class(fit_d)=="try-error"){return(ori_p);}
		pv <- pgpd(ob-t_score,mu=0,beta=fit_d@fit$par.ests["beta"],xi=fit_d@fit$par.ests["xi"]);
		adj_p <- ori_p*(1-pv[1]);
		return(as.numeric(adj_p));
		}'''
		pareto=r(Rcode)
###############KEGG
		pathway_name=list()
		p_value_vector=list()
		score_vector=list()
		node_number=list()
		databasefile=dir+"/data/"+ "KEGG_pathway_index.txt"
		for line in open(databasefile):
			line=line.strip('\n').split('\t')
			pathway_name.append(line[0])
			score_total=np.mean(all_score[list(map(int, line[1:]))],0)
			score_vector.append(score_total[0])
			p_value=len(score_total[score_total>score_total[0]])/random_number
			if p_value<=p_threshold:
				args[0]=score_total[0]
				p_value=np.array(pareto(robjects.FloatVector(score_total[1:(random_number+1)]),robjects.FloatVector(args)))[0]

			p_value_vector.append(p_value)
			node_number.append(len(line)-1)

		kegg=pd.DataFrame({"pathway_name":pathway_name,"PS":score_vector,"p_value":p_value_vector,"node_number":node_number})
		kegg['Source']="KEGG"
		websitelink=dir+"/data/KEGG_website.txt"
		link=pd.read_csv(websitelink,sep="\t")
		kegg=pd.merge(kegg,link,how='left',on='pathway_name')
		kegg=kegg.fillna("null")
###############Reactome
		pathway_name=list()
		p_value_vector=list()
		score_vector=list()
		node_number=list()
		databasefile=dir+"/data/"+ "Reactome_pathway_index.txt"
		for line in open(databasefile):
			line=line.strip('\n').split('\t')
			pathway_name.append(line[0])
			score_total=np.mean(all_score[list(map(int, line[1:]))],0)
			score_vector.append(score_total[0])
			p_value=len(score_total[score_total>score_total[0]])/random_number
			if p_value<=p_threshold:
				args[0]=score_total[0]
				p_value=np.array(pareto(robjects.FloatVector(score_total[1:(random_number+1)]),robjects.FloatVector(args)))[0]

			p_value_vector.append(p_value)
			node_number.append(len(line)-1)

		Reactome=pd.DataFrame({"pathway_name":pathway_name,"PS":score_vector,"p_value":p_value_vector,"node_number":node_number})
		Reactome['Source']="Reactome"
		websitelink=dir+"/data/Reactome_website.txt"
		link=pd.read_csv(websitelink,sep="\t")
		Reactome=pd.merge(Reactome,link,how='left',on='pathway_name')
		Reactome=Reactome.fillna("null")
###########Pid
		pathway_name=list()
		p_value_vector=list()
		score_vector=list()
		node_number=list()
		databasefile=dir+"/data/"+ "PID_pathway_index.txt"
		for line in open(databasefile):
			line=line.strip('\n').split('\t')
			pathway_name.append(line[0])
			score_total=np.mean(all_score[list(map(int, line[1:]))],0)
			score_vector.append(score_total[0])
			p_value=len(score_total[score_total>score_total[0]])/random_number
			if p_value<=p_threshold:
				args[0]=score_total[0]
				p_value=np.array(pareto(robjects.FloatVector(score_total[1:(random_number+1)]),robjects.FloatVector(args)))[0]

			p_value_vector.append(p_value)
			node_number.append(len(line)-1)

		Pid=pd.DataFrame({"pathway_name":pathway_name,"PS":score_vector,"p_value":p_value_vector,"node_number":node_number})
		Pid['Source']="PID"
		websitelink=dir+"/data/PID_website.txt"
		link=pd.read_csv(websitelink,sep="\t")
		Pid=pd.merge(Pid,link,how='left',on='pathway_name')
		Pid=Pid.fillna("null")
###########Netpath
		pathway_name=list()
		p_value_vector=list()
		score_vector=list()
		node_number=list()
		databasefile=dir+"/data/"+ "NetPath_pathway_index.txt"
		for line in open(databasefile):
			line=line.strip('\n').split('\t')
			pathway_name.append(line[0])
			score_total=np.mean(all_score[list(map(int, line[1:]))],0)
			score_vector.append(score_total[0])
			p_value=len(score_total[score_total>score_total[0]])/random_number
			if p_value<=p_threshold:
				args[0]=score_total[0]
				p_value=np.array(pareto(robjects.FloatVector(score_total[1:(random_number+1)]),robjects.FloatVector(args)))[0]

			p_value_vector.append(p_value)
			node_number.append(len(line)-1)

		Netpath=pd.DataFrame({"pathway_name":pathway_name,"PS":score_vector,"p_value":p_value_vector,"node_number":node_number})
		Netpath['Source']="NetPath"
		websitelink=dir+"/data/NetPath_website.txt"
		link=pd.read_csv(websitelink,sep="\t")
		Netpath=pd.merge(Netpath,link,how='left',on='pathway_name')
		Netpath=Netpath.fillna("null")
#########Inoh
		pathway_name=list()
		p_value_vector=list()
		score_vector=list()
		node_number=list()
		databasefile=dir+"/data/"+ "INOH_pathway_index.txt"
		for line in open(databasefile):
			line=line.strip('\n').split('\t')
			pathway_name.append(line[0])
			score_total=np.mean(all_score[list(map(int, line[1:]))],0)
			score_vector.append(score_total[0])
			p_value=len(score_total[score_total>score_total[0]])/random_number
			if p_value<=p_threshold:
				args[0]=score_total[0]
				p_value=np.array(pareto(robjects.FloatVector(score_total[1:(random_number+1)]),robjects.FloatVector(args)))[0]

			p_value_vector.append(p_value)
			node_number.append(len(line)-1)

		Inoh=pd.DataFrame({"pathway_name":pathway_name,"PS":score_vector,"p_value":p_value_vector,"node_number":node_number})
		Inoh['Source']="INOH"
		websitelink=dir+"/data/INOH_website.txt"
		link=pd.read_csv(websitelink,sep="\t")
		Inoh=pd.merge(Inoh,link,how='left',on='pathway_name')
		Inoh=Inoh.fillna("null")
#########Humancyc
		pathway_name=list()
		p_value_vector=list()
		score_vector=list()
		node_number=list()
		databasefile=dir+"/data/"+ "HumanCyc_pathway_index.txt"
		for line in open(databasefile):
			line=line.strip('\n').split('\t')
			pathway_name.append(line[0])
			score_total=np.mean(all_score[list(map(int, line[1:]))],0)
			score_vector.append(score_total[0])
			p_value=len(score_total[score_total>score_total[0]])/random_number
			if p_value<=p_threshold:
				args[0]=score_total[0]
				p_value=np.array(pareto(robjects.FloatVector(score_total[1:(random_number+1)]),robjects.FloatVector(args)))[0]

			p_value_vector.append(p_value)
			node_number.append(len(line)-1)

		Humancyc=pd.DataFrame({"pathway_name":pathway_name,"PS":score_vector,"p_value":p_value_vector,"node_number":node_number})
		Humancyc['Source']="HumanCyc"
		websitelink=dir+"/data/HumanCyc_website.txt"
		link=pd.read_csv(websitelink,sep="\t")
		Humancyc=pd.merge(Humancyc,link,how='left',on='pathway_name')
		Humancyc=Humancyc.fillna("null")
##########Panther
		pathway_name=list()
		p_value_vector=list()
		score_vector=list()
		node_number=list()
		databasefile=dir+"/data/"+ "PANTHER_pathway_index.txt"
		for line in open(databasefile):
			line=line.strip('\n').split('\t')
			pathway_name.append(line[0])
			score_total=np.mean(all_score[list(map(int, line[1:]))],0)
			score_vector.append(score_total[0])
			p_value=len(score_total[score_total>score_total[0]])/random_number
			if p_value<=p_threshold:
				args[0]=score_total[0]
				p_value=np.array(pareto(robjects.FloatVector(score_total[1:(random_number+1)]),robjects.FloatVector(args)))[0]

			p_value_vector.append(p_value)
			node_number.append(len(line)-1)

		Panther=pd.DataFrame({"pathway_name":pathway_name,"PS":score_vector,"p_value":p_value_vector,"node_number":node_number})
		Panther['Source']="PANTHER"
		websitelink=dir+"/data/PANTHER_website.txt"
		link=pd.read_csv(websitelink,sep="\t")
		Panther=pd.merge(Panther,link,how='left',on='pathway_name')
		Panther=Panther.fillna("null")
##########Wikipathways
		pathway_name=list()
		p_value_vector=list()
		score_vector=list()
		node_number=list()
		databasefile=dir+"/data/"+ "WikiPathways_pathway_index.txt"
		for line in open(databasefile):
			line=line.strip('\n').split('\t')
			pathway_name.append(line[0])
			score_total=np.mean(all_score[list(map(int, line[1:]))],0)
			score_vector.append(score_total[0])
			p_value=len(score_total[score_total>score_total[0]])/random_number
			if p_value<=p_threshold:
				args[0]=score_total[0]
				p_value=np.array(pareto(robjects.FloatVector(score_total[1:(random_number+1)]),robjects.FloatVector(args)))[0]

			p_value_vector.append(p_value)
			node_number.append(len(line)-1)

		Wikipathways=pd.DataFrame({"pathway_name":pathway_name,"PS":score_vector,"p_value":p_value_vector,"node_number":node_number})
		Wikipathways['Source']="WikiPathways"
		websitelink=dir+"/data/WikiPathways_website.txt"
		link=pd.read_csv(websitelink,sep="\t")
		Wikipathways=pd.merge(Wikipathways,link,how='left',on='pathway_name')
		Wikipathways=Wikipathways.fillna("null")
#############
		out_result=pd.concat([kegg,Inoh,Reactome,Netpath,Pid,Panther,Wikipathways,Humancyc])
		out_result=out_result.sort_values(['p_value','PS'],ascending=[True,False])
		adj_p=adjust_pvalues(out_result.p_value)
		out_result['FDR']=adj_p
		order=['pathway_name','Source','PS','p_value','FDR','node_number',"website"]
		out_result=out_result[order]
		out_result['p_value']=out_result['p_value'].map(lambda x:('%.3e')%x)
		out_result['FDR']=out_result['FDR'].map(lambda x:('%.3e')%x)
		out_result['PS']=out_result['PS'].map(lambda x:('%.4f')%x)
		out_result=out_result[order]
		out_result.to_csv(outfile,index=False,header=True,sep="\t")
		if sig_method=="FDR":
			used= out_result.loc[out_result['FDR'].astype('float')< float(sig_value)]
		else:
			used= out_result.loc[out_result['p_value'].astype('float')< float(sig_value)]

		used.to_csv(sigfile,index=False,header=True,sep="\t")
		path_svg= path
		out_svg_1 = path_svg + "/bar_20.svg"
		out_png_1 = path + "/bar.pdf"
		out_bubble_svg_1 = path_svg + "/bubble_20.svg"
		out_bubble_png_1 = path + "/bubble.pdf"
		plt_bubble_bar_svg(outfile,out_bubble_svg_1,out_svg_1,20,"All")
		cairosvg.svg2pdf(url=out_bubble_svg_1,write_to=out_bubble_png_1)
		cairosvg.svg2pdf(url=out_svg_1,write_to=out_png_1)
		os.remove(out_bubble_svg_1)
		os.remove(out_svg_1)
		index=0.3
		index_list=[0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9]
		for i in index_list:
			cytoscape(i,sig_method,sig_value,"All",dir)
			cytoscape_file=Path(path_svg + "/network_temp_"+str(i)+".null")
			if cytoscape_file.is_file():
				os.remove(cytoscape_file)
				break


		email_content ='Congratulations! Your CTpathway job has been completed. '
		print(email_content)
	except Exception as e:
		error_log=path + '/error.txt'
		traceback.print_exc(file=open(error_log,'a+'))
		email_content = 'Sorry, your CTpathway job is failed! Please check your input file!'
		print(email_content)
		#shutil.rmtree(path)

if __name__=='__main__':
	import os
	import argparse
	parser = argparse.ArgumentParser()
	parser.add_argument('file', help='Gene differential expression profile, which should include gene id, one or both of log2FC (Fold Change) value and P value. Notice: input file without head line is commanded.')
	parser.add_argument('gene_id', help='Gene id type corresponds to gene differential expression profile. Available options are "symbol" or "entrez_id".')
	parser.add_argument('data_type', help='Data type corresponds to gene differential expression profile. Available options are "Log2_FC", "P_value" or "Both".')
	parser.add_argument('--pathway_database', default='KEGG', help='Optional. Select pathway database used in this analysis.  Availale options are "KEGG", "Reactome", "PANTHER", "HumanCyc", "INOH", "NetPath", "PID", "WikiPathways" or "All". [default= "KEGG"].')
	parser.add_argument('--permutation_number', default='1000', help='Optional. Permutation number is used to identify significant pathways. [default= 1000]. Notice: the larger permutation number is, the longer time consumes.')
	parser.add_argument('--r_value', default='0.7', help='Optional. Restart coefficient corresponds to multi-RWR (multiple random walk with restart) algorithm. Available options are "0.1", "0.2", "0.3", "0.4", "0.5", "0.6", "0.7", "0.8" or "0.9". [default= 0.7]. Notice: the larger value represents more dependency on differential expression and less dependency on GPCM, which is based on prior biological knowledge.')
	parser.add_argument('--significant_method', default='FDR', help='Optional. Significant method is to choose P value or FDR to screen significant pathways. Availale options are "FDR" or "P_value". [default= "FDR"].')
	parser.add_argument('--significant_threshold', default='0.01', help='Optional. Significant threshold is threshold of significant pathways. [default= 0.01].')
	args = parser.parse_args()
	print(args.file)
	print(args.gene_id)
	print(args.data_type)
	print(args.significant_method)
	print(args.significant_threshold)
	print(args.r_value)
	print(args.pathway_database)
	print(args.permutation_number)
	dir=os.path.split(os.path.realpath(__file__))[0]
	outdir=dir+"/output"
	if os.path.exists(outdir):
		shutil.rmtree(outdir)
		os.mkdir(outdir)
	else:
		os.mkdir(outdir)

	rwrfile=dir+'/data/sparse_filter.3.all_final_'+args.r_value+'.txt.npz'
	if args.gene_id == "symbol":
		if args.pathway_database != "All": 
			calculate_pathway_score_symbol(args.file,args.pathway_database,rwrfile,args.permutation_number,args.data_type,args.significant_method,args.significant_threshold,dir)
		else:
			calculate_pathway_score_symbol_all(args.file,args.pathway_database,rwrfile,args.permutation_number,args.data_type,args.significant_method,args.significant_threshold,dir)
	
	if args.gene_id == "entrez_id":
		if args.pathway_database != "All": 
			calculate_pathway_score(args.file,args.pathway_database,rwrfile,args.permutation_number,args.data_type,args.significant_method,args.significant_threshold,dir)
		else:
			calculate_pathway_score_all(args.file,args.pathway_database,rwrfile,args.permutation_number,args.data_type,args.significant_method,args.significant_threshold,dir)