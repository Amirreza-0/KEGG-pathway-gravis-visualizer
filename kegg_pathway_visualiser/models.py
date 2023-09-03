import os

from django.core.files import File
from django.db import models

# Create your models here.
class KEGGML_Model(models.Model):

    def name1(instance, filename):
        return os.path.join('kegg_pathway_visualiser/KEGGML_files/', filename)

    name = models.CharField(max_length=50)

    file = models.FileField()
    uploaded_at = models.DateTimeField(auto_now_add=True)




