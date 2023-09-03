from django import forms
from .models import KEGGML_Model
# Create your forms here.

class KEGGML_form(forms.ModelForm):
    # list of target genes
    target_genes = forms.CharField(max_length=1000)

    target_genes.widget.attrs.update({'class': 'form-control', 'placeholder': 'TP53,MDM2'})
    # not required
    target_genes.required = False
    target_relations = forms.CharField(max_length=1000)
    target_relations.widget.attrs.update({'class': 'form-control', 'placeholder': 'activation,inhibition'})
    # not required
    target_relations.required = False
    target_generations = forms.IntegerField()
    # not required
    target_generations.required = False

    class Meta:
        model = KEGGML_Model
        fields = ('name', 'file', 'target_genes', 'target_relations', 'target_generations')

    # converting the target genes and relations to list
    def clean_target_genes(self):
        target_genes = self.cleaned_data['target_genes']
        if ("," in target_genes):
            target_genes = target_genes.split(',')
        else:
            target_genes = [target_genes]
        return target_genes

    def clean_target_relations(self):
        target_relations = self.cleaned_data['target_relations']
        if ("," in target_relations):
            target_relations = target_relations.split(',')
        else:
            target_relations = [target_relations]
        return target_relations


