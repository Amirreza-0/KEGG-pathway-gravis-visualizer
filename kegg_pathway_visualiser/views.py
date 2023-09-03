from django.shortcuts import render, redirect
from django.http import HttpResponse
from .forms import KEGGML_form
from .models import KEGGML_Model
from .multiprocessing_symbols_extraction_ import extract_symbols
from .Build_Adjacency_Matrix import build_adjacency_matrix
from .plotly_html_generator import generate_plotly_html
from .Gravis import generate_gravis_html


def process(target_genes, target_relations, target_generations):
    file_name = KEGGML_Model.objects.last().file.name
    symbols_path = extract_symbols(file_name)
    matrix_path = build_adjacency_matrix(file_name, symbols_path=symbols_path)
    html, img = generate_gravis_html(matrix_path, symbols_path=symbols_path, target_genes=target_genes,
                                target_relations=target_relations, target_generations=target_generations)

    return html, img
# Create your views here.
def index(request):
    # loading the previous input file from KEGGML model
    #file_names = KEGGML_Model.objects.all()


    # saving the input file in KEGML model and redirecting to the loading page
    if request.method == 'POST':
        form = KEGGML_form(request.POST, request.FILES)
        if form.is_valid():
            form.save()
            target_genes = form.clean_target_genes()
            target_relations = form.clean_target_relations()
            target_generations = form.cleaned_data['target_generations']
            html, img = process(target_genes, target_relations, target_generations)
            return render(request, 'loading_results.html', context={"html": html, "img": img})

        if form.errors:
            print(form.errors)
            return redirect('index')

    # if a GET (or any other method) we'll create a blank form
    else:
        form = KEGGML_form()

    return render(request, 'index.html', context={"form": form})




"""
def loading_results(request):
    # loading the input file from KEGGML model and extracting symbols
    if request.method == 'GET':
        file_name = KEGGML_Model.objects.last().file.name
        symbols_path = extract_symbols(file_name)
        matrix_path = build_adjacency_matrix(file_name, symbols_path=symbols_path)
        html = generate_gravis_html(matrix_path, symbols_path=symbols_path, target_genes=target_genes,
                                    target_relations=target_relations, target_generations=target_generations)
        return render(request, 'loading_results.html', context={"html": html})

    return render(request, 'loading_results.html')
"""