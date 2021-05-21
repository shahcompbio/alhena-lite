from flask import Flask, session, jsonify
from flask_session import Session
import redis
import pandas as pd
import numpy as np
from scgenome.loaders.qc import load_qc_data

app = Flask(__name__)
app.config['SESSION_TYPE'] = 'redis'
app.config['SESSION_PERMANENT'] = False
app.config['SESSION_USE_SIGNER'] = True
app.config['SESSION_REDIS'] = redis.from_url('redis://localhost:6379')
server_session = Session(app)

app.secret_key = "alhenalite"


@app.route('/')
def hello_world():
    return 'Hello World!'


@app.route('/api/<path:directory>/')
def api(directory=None):
    # given directory, return cell, seg, gc data
    session.clear()
    if directory:
        data = load_qc_data("/" + directory)
        session['bins'] = get_bins_data(data)
        gc_bias = get_gc_bias_data(data)
        segs = get_segs_data(data)
        qc = get_qc_data(data)

        grouped_segs = segs.groupby(['id'])

        qc_records = []
        for record in qc.to_dict(orient="records"):
            clean_nans(record)
            cell_segs = []
            for seg_record in grouped_segs.get_group(record['id']).to_dict(orient='record'):
                clean_nans(seg_record)
                cell_segs.append(seg_record)

            record['segs'] = cell_segs
            qc_records.append(record)

        return jsonify({
            'gc_bias': gc_bias.to_dict(orient='record'),
            'cells': qc_records
        })


@app.route('/cell_ids')
def cell_ids():
    return jsonify(session['qc']['id'].tolist())


@app.route('/bin/<cell_id>')
def bin(cell_id):
    bins = session['bins']
    bins = bins.loc[bins['id'] == cell_id]

    records = []
    for record in bins.to_dict(orient="records"):
        clean_nans(record)
        records.append(record)
    return jsonify(records)


def get_qc_data(hmmcopy_data):
    data = hmmcopy_data['annotation_metrics']
    data['percent_unmapped_reads'] = data["unmapped_reads"] / data["total_reads"]
    data['is_contaminated'] = data['is_contaminated'].apply(
        lambda a: {True: 'true', False: 'false'}[a])
    data = data.rename(columns={'cell_id': 'id'})
    return data


def get_segs_data(hmmcopy_data):
    data = hmmcopy_data['hmmcopy_segs'].copy()
    data['chrom_number'] = create_chrom_number(data['chr'])
    data = data.rename(columns={'cell_id': 'id'})
    return data


def get_bins_data(hmmcopy_data):
    data = hmmcopy_data['hmmcopy_reads'].copy()
    data['chrom_number'] = create_chrom_number(data['chr'])
    data = data.rename(columns={'cell_id': 'id'})

    return data


def get_gc_bias_data(hmmcopy_data):
    data = hmmcopy_data['gc_metrics']

    gc_cols = list(range(101))
    gc_bias_df = pd.DataFrame(columns=['cell_id', 'gc_percent', 'value'])
    for n in gc_cols:
        new_df = data.loc[:, ['cell_id', str(n)]]
        new_df.columns = ['cell_id', 'value']
        new_df['gc_percent'] = n
        gc_bias_df = gc_bias_df.append(new_df, ignore_index=True)

    gc_bias_df = gc_bias_df.rename(columns={'cell_id': 'id'})

    return gc_bias_df


chr_prefixed = {str(a): '0' + str(a) for a in range(1, 10)}


def create_chrom_number(chromosomes):
    chrom_number = chromosomes.map(lambda a: chr_prefixed.get(a, a))
    return chrom_number


def clean_nans(record):
    """Delete any fields that contain nan values"""
    floats = [field for field in record if isinstance(record[field], float)]
    for field in floats:
        if np.isnan(record[field]):
            del record[field]
