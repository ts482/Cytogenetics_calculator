from flask import Flask, request
from flask_cors import CORS
import Extractor as cyto

app = Flask(__name__)
CORS(app)

#setup
abnormalities = cyto.base_extraction()
props = cyto.setup(abnormalities)



@app.route("/")
def hello_world():
    return "<p>Hello, World!</p>"

@app.route('/karyotype/extract', methods=['POST'])
def karyotype_api():
    input = request.json
    if 'fish' not in input:
        input['fish'] = None
    extracted = cyto.extract_from_string(input['karyotype_string'], props, fish = input['fish'])
    return extracted

# TODO bulk API to extract multiple strings at once

# could easily add an endpoint that also allows the abnormalities for extraction to be specified
# rather than using the app-level default 

@app.route('/karyotype/ping', methods=['POST'])
def ping():
    input = request.json
    return input