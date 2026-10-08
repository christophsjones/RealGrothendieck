"""Independent traversal of the delivered finite-domain proof trees.

No numerical optimization, floating-point upper bound, or discovery stopping
condition is used for acceptance. Every flow is regenerated and its feasible
flow / equal-capacity cut witness is checked by the integer C++ core.
"""
if not __debug__:
    raise RuntimeError('Run certificate verification without Python -O')
import os
os.environ['OPENBLAS_NUM_THREADS']='1'
os.environ['OMP_NUM_THREADS']='1'
from pathlib import Path
from fractions import Fraction as F
import json,hashlib,time
import numpy as np
from graph import build,roof,cap,symmetry,validate
S=1<<44
D=65536


def unique_object(pairs):
    out={}
    for k,v in pairs:
        if k in out: raise ValueError('Duplicate JSON key: '+k)
        out[k]=v
    return out


def load_json(path):
    return json.loads(Path(path).read_text(),object_pairs_hook=unique_object)


def sha(path): return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def require(test,message):
    if not test: raise ValueError(message)


def exact_int(x):
    require(type(x) is int,'Expected an exact JSON integer')
    return x


def eliminate(B,a,c,indices,signs):
    """Substitute fixed signs, with a checked global int64 safety budget."""
    validate(B,a,c)
    indices=np.asarray(indices,dtype=np.int64)
    signs=np.asarray(signs,dtype=np.int64)
    require(len(indices)==len(signs) and len(set(map(int,indices)))==len(indices),'Invalid fixed coordinates')
    require(bool(np.all((indices>=0)&(indices<len(a)))),'Fixed index out of range')
    require(bool(np.all(np.abs(signs)==1)),'Fixed values are not signs')
    mask=np.ones(len(a),dtype=bool);mask[indices]=False
    rest=np.flatnonzero(mask)
    twice=int(signs@B[np.ix_(indices,indices)]@signs)
    require(twice%2==0,'Odd doubled diagonal contribution')
    cc=int(c)+int(a[indices]@signs)+twice//2
    aa=np.ascontiguousarray(a[rest]+B[np.ix_(rest,indices)]@signs)
    BB=np.ascontiguousarray(B[np.ix_(rest,rest)])
    validate(BB,aa,cc)
    return BB,aa,cc,rest


def field_coordinate(a):
    den=sum(abs(int(x)) for x in a)
    require(den>0,'Zero coordinate normalizer')
    result=[]
    for z in a:
        z=int(z);n=abs(z)*(1<<24)
        result.append((1 if z>=0 else -1)*((2*n+den)//(2*den)))
    return np.ascontiguousarray(result,dtype=np.int64)


def checked_node_map(nodes):
    require(type(nodes) is list and len(nodes)>0,'Missing proof nodes')
    table={}
    for node in nodes:
        require(type(node) is dict,'Invalid node')
        k=exact_int(node['id']);require(k>=0 and k not in table,'Duplicate or negative node id')
        table[k]=node
    require(0 in table,'Missing root node')
    return table


def replay_binary(B,a,c,target,record,allow_prefix=False,mean_coordinate=None):
    require(record['status']=='EXACT_BINARY_UPPER_CERTIFIED','Incomplete binary proof')
    require(record['active_count']==0,'Active binary nodes remain')
    require(record['target_integer']==target,'Binary target mismatch')
    require(record['all_flows_checked'] is True,'Missing flow check flag')
    require(record['node_count']==len(record['nodes']),'Binary node count mismatch')
    validate(B,a,c)
    prefix=record['reflection_prefix']
    if prefix:
        require(allow_prefix and mean_coordinate is not None,'Unjustified binary symmetry prefix')
        N=len(a)//2
        require(len(a)==2*N and N%2==1 and prefix==[[N//2,1]],'Invalid reflection prefix')
        require(symmetry(B,a,mean_coordinate),'Reflection symmetry not satisfied')
        B,a,c,rest=eliminate(B,a,c,[N//2],[1])
        ids=np.arange(2*N,dtype=np.int32)[rest]
    else:
        ids=np.arange(len(a),dtype=np.int32)
    table=checked_node_map(record['nodes']);seen=set()
    stack=[(0,B,a,c,ids)];maximum=None;flow_count=0;leaves=0
    while stack:
        nodeid,B,a,c,ids=stack.pop()
        require(nodeid in table and nodeid not in seen,'Missing/repeated binary node')
        seen.add(nodeid);node=table[nodeid]
        require(node['dimension']==len(a),'Binary dimension mismatch')
        if len(a): bound,labels,flow=roof(B,a,c)
        else: bound,labels,flow=c,np.zeros(0,dtype=np.int32),0
        flow_count+=int(len(a)>0)
        require(node['bound']==bound and node['flow']==flow,'Binary flow arithmetic mismatch')
        if node['kind']=='leaf':
            require(bound<=target,'Unproved binary terminal bound')
            maximum=bound if maximum is None else max(maximum,bound);leaves+=1
            require('children' not in node,'Terminal binary node has children')
            continue
        require(node['kind']=='split','Unaccepted binary node kind')
        fix=np.flatnonzero(labels!=0)
        expected=[[int(ids[j]),int(labels[j])] for j in fix]
        require(node['persistent']==expected,'Unjustified persistent assignments')
        BB,aa,cc,rest=eliminate(B,a,c,fix,labels[fix]);newids=ids[rest]
        v=exact_int(node['variable']);where=np.flatnonzero(newids==v)
        require(len(where)==1,'Binary split variable not free')
        j=int(where[0]);children=node['children']
        require(children==[2*nodeid+1,2*nodeid+2],'Invalid binary child identifiers')
        for sign,child in [(-1,children[0]),(1,children[1])]:
            BBB,aaa,ccc,rr=eliminate(BB,aa,cc,[j],[sign])
            stack.append((child,BBB,aaa,ccc,newids[rr]))
    require(seen==set(table),'Extraneous binary proof nodes')
    require(maximum is not None and maximum==record['certified_upper_integer'],'Binary assembly mismatch')
    return maximum,dict(integer_flows=flow_count,binary_nodes=len(table),binary_leaves=leaves)


def regenerated_caps(B,a,c,vs,box,prices,prefix):
    require(type(prices) is list and len(prices)==2,'Missing cap prices')
    error=0
    for k,p in enumerate(prices):
        exact_int(p);require(0<=p<=1<<30,'Invalid/nonpositive cap price')
        if p:
            B,a,c,e=cap(B,a,c,vs[k],p,box[0][k],box[1][k]);error+=e
    B,a,c,rest=eliminate(B,a,c,[prefix],[1])
    return B,a,c,error


def replay_boxes(B,a,c,mean,target,record):
    require(record['status']=='UNIVERSAL_FIBRE_CERTIFIED','Incomplete rectangle proof')
    require(record['pending']==[],'Pending rectangles remain')
    require(record['target_integer']==target,'Rectangle target mismatch')
    N=len(a)//2
    require(len(a)==2*N and N%2==1,'Invalid grid size for reflection prefix')
    prefix=N//2
    require(record['prefix']==[prefix,1],'Invalid external reflection prefix')
    require(symmetry(B,a,mean),'Original graph lacks reflection symmetry')
    vs=[field_coordinate(a),mean]
    norms=[sum(abs(int(z)) for z in v) for v in vs]
    require(all(0<V<=1<<30 for V in norms),'Coordinate norm out of range')
    require(record['coordinates']==[v.tolist() for v in vs],'Coordinate vectors changed')
    require(record['coordinate_normalizers']==norms,'Coordinate normalizers changed')
    table=checked_node_map(record['nodes']);seen=set();maximum=None
    stack=[(0,[-D,-D],[D,D])]
    stats=dict(rectangles=len(table),empty_rectangles=0,flow_rectangles=0,branch_rectangles=0,integer_flows=0,binary_nodes=0,binary_leaves=0)
    while stack:
        nodeid,lo,hi=stack.pop()
        require(nodeid in table and nodeid not in seen,'Missing/repeated rectangle node')
        seen.add(nodeid);node=table[nodeid]
        require(node['lo']==lo and node['hi']==hi,'Rectangle gap or endpoint mismatch')
        require(all(type(z) is int for z in lo+hi),'Nonintegral rectangle endpoint')
        require(all(-D<=lo[k]<hi[k]<=D for k in range(2)),'Invalid rectangle')
        kind=node['kind']
        if kind=='split':
            k=exact_int(node['coordinate']);m=exact_int(node['at'])
            require(k in [0,1] and lo[k]<m<hi[k],'Invalid rectangle split')
            children=node['children'];require(children==[2*nodeid+1,2*nodeid+2],'Invalid rectangle children')
            h0=hi.copy();h0[k]=m;l1=lo.copy();l1[k]=m
            stack.extend([(children[0],lo.copy(),h0),(children[1],l1,hi.copy())]);continue
        require('children' not in node,'Terminal rectangle has children')
        if kind=='empty':
            facet=node['facet'];require(type(facet) is list and len(facet)==3,'Missing facet witness')
            alpha,beta,h=map(exact_int,facet);require(alpha!=0 or beta!=0,'Zero facet')
            support=alpha*int(vs[0][prefix])+beta*int(vs[1][prefix])
            support+=sum(abs(alpha*int(vs[0][j])+beta*int(vs[1][j])) for j in range(len(a)) if j!=prefix)
            require(h==support,'Invalid exact coordinate support')
            minimum=alpha*norms[0]*(lo[0] if alpha>=0 else hi[0])+beta*norms[1]*(lo[1] if beta>=0 else hi[1])
            require(minimum>h*D,'Rectangle is not excluded by its support inequality')
            stats['empty_rectangles']+=1;continue
        require(kind in ['flow','branch'],'Unaccepted rectangle leaf')
        cert=node['certificate'];BB,aa,cc,e=regenerated_caps(B,a,c,vs,(lo,hi),cert['prices'],prefix)
        require(cert['rounding_error_units']==e,'Cap rounding error mismatch')
        if kind=='flow':
            bd,labels,flow=roof(BB,aa,cc)
            require(cert['bound']==bd and cert['flow']==flow,'Rectangle flow arithmetic mismatch')
            stats['integer_flows']+=1;stats['flow_rectangles']+=1
        else:
            bd,extra=replay_binary(BB,aa,cc,target-e,cert['tree'])
            for k,v in extra.items():stats[k]+=v
            stats['branch_rectangles']+=1
        upper=bd+e
        require(cert['upper']==upper and upper<=target,'Unproved rectangle upper bound')
        maximum=upper if maximum is None else max(maximum,upper)
    require(seen==set(table),'Extraneous rectangle nodes')
    require(maximum is not None and maximum==record['max_upper_integer'],'Rectangle assembly mismatch')
    return maximum,stats


def regenerate_analytic(root,row,out):
    from analytic import Common,read_row,check_pd,curvature,grid_nodes,integral_error,hat_inputs,save_input,upper
    root=Path(root);out=Path(out);out.mkdir(parents=True,exist_ok=True)
    source=root/'baseline/lower_push_20260906/multiple_caps/best_candidate/frozen_rows.json'
    meta=load_json(root/'inputs'/f'row{row}_analytic.json')
    require(meta['row']==row and meta['source_sha256']==sha(source),'Analytic source linkage mismatch')
    require(meta['precision_bits']==384,'Unexpected analytic precision')
    require(meta.get('method')=='quadratic_bernstein','Wrong interpolation method')
    common=Common();M,ell=read_row(source,row);pd=check_pd(M)
    require(F(pd)>0 and F(meta['pd_min_pivot_lower'])>0,'Positive-definiteness failed')
    vals,tail=curvature(common,M,ell)
    require(all(upper(v)<=F(q) for v,q in zip(vals,meta['curvature_bounds'])) and len(meta['curvature_bounds'])==36,'Curvature upper bounds are insufficient')
    require(upper(tail)<=F(meta['tail_error']),'Tail error is insufficient')
    counts=meta['counts'];nodes=grid_nodes(counts)
    require(len(nodes)==meta['n'],'Grid size mismatch')
    e1=integral_error(common,counts,vals,tail)
    require(e1<=F(meta['interpolation_and_tail_error']),'Interpolation error is insufficient')
    Bint,mint,eint,e2,coords,mass=hat_inputs(common,nodes,M,ell)
    require(e2<=F(meta['fixed_point_input_error']),'Fixed-point input error is insufficient')
    require(len(coords)==len(meta['coordinate_errors']) and all(F(q)<=F(r) for q,r in zip(coords,meta['coordinate_errors'])),'Coordinate error mismatch')
    fresh=out/f'row{row}.bin';hash0=save_input(fresh,Bint,mint,eint)
    original=root/'inputs'/f'row{row}.bin'
    require(hash0==meta['input_sha256']==sha(original),'Fixed-point input regeneration mismatch')
    require(fresh.read_bytes()==original.read_bytes(),'Fixed-point input byte mismatch')
    return meta,fresh


def verify_row(root,row,out=None,regenerate=True):
    tic=time.time();root=Path(root);out=Path(out or root/'replay');out.mkdir(parents=True,exist_ok=True)
    if regenerate:meta,inp=regenerate_analytic(root,row,out/'inputs')
    else:
        # Used only by isolated malformed-evidence tests, never final acceptance.
        meta=load_json(root/'inputs'/f'row{row}_analytic.json');inp=root/'inputs'/f'row{row}.bin'
        require(sha(inp)==meta['input_sha256'],'Input hash mismatch')
    B,a,c,qe,mean=build(inp)
    targets=load_json(root/'inputs/target_allowances.json')['rows']
    require(len(targets)==35 and targets[row]['index']==row,'Target row mismatch')
    allow=F(targets[row]['allowed_fibre_upper'])
    error=F(meta['interpolation_and_tail_error'])+F(meta['fixed_point_input_error'])+qe
    scaled=(allow-error)*S;target=scaled.numerator//scaled.denominator
    p=root/'certificates'/f'row{row}.json';data=load_json(p)
    require(data['row']==row and data['input_sha256']==meta['input_sha256'],'Finite certificate input mismatch')
    if data['method']=='binary':
        maximum,stats=replay_binary(B,a,c,target,data['proof'],allow_prefix=True,mean_coordinate=mean)
    elif data['method']=='boxes':
        record=data['proof'];require(record['row']==row,'Rectangle row mismatch')
        require(record['input_sha256']==meta['input_sha256'],'Rectangle input mismatch')
        require(F(record['analytic_and_input_error'])==error and F(record['allowed'])==allow,'Rectangle rational assembly mismatch')
        maximum,stats=replay_boxes(B,a,c,mean,target,record)
        require(F(record['certified_fibre_upper'])==F(maximum,S)+error,'Stored fibre upper mismatch')
    else:raise ValueError('Unknown finite certificate method')
    upper_bound=F(maximum,S)+error
    require(upper_bound<=allow,'Universal fibre allowance exceeded')
    result=dict(status='UNIVERSAL_FIBRE_REPLAY_PASSED' if regenerate else 'FINITE_GRAPH_REPLAY_ONLY',row=row,method=data['method'],source_sha256=meta['source_sha256'],input_sha256=meta['input_sha256'],certificate_sha256=sha(p),analytic_regenerated=regenerate,precision_bits=384,grid_nodes_per_channel=len(a)//2,interpolation_and_tail_error=meta['interpolation_and_tail_error'],fixed_point_input_error=meta['fixed_point_input_error'],graph_quantization_error=str(qe),total_error=str(error),finite_upper_integer=maximum,finite_scale=S,certified_fibre_upper=str(upper_bound),allowed=str(allow),strict_slack=str(allow-upper_bound),statistics=stats,seconds=time.time()-tic)
    (out/f'row{row}.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:result[k] for k in ['status','row','method','statistics','seconds']}),flush=True)
    return result

if __name__=='__main__':
    import argparse
    ap=argparse.ArgumentParser();ap.add_argument('--root',default=str(Path(__file__).resolve().parents[1]));ap.add_argument('--row',type=int,required=True);ap.add_argument('--out');args=ap.parse_args()
    verify_row(args.root,args.row,args.out)
