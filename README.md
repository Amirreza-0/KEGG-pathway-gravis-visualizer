# KEGG-KGML pathway visualizer
 
### how to run the django project
1. create a virtual environment using `python -m venv venv` and activate it using `venv\Scripts\activate`
2. install the requirements.txt
3. make migrations using `python manage.py makemigrations`
4. migrate to build the initial database using `python manage.py migrate`
5. run the django project using `python manage.py runserver`

then you can start using the project by following the given link in the terminal by default it is `http://127.0.0.1:8000/

The test data that has been used for generating results can be found on KEGG website by following this link:
https://www.kegg.jp/pathway/hsa05200

any other standard KGML files from KEGG website should work.

to generate exact same results shown below you can give these parameters targets_genes: TP53,MDM2 target relations: activation,inhibition,expression,repression,indirect effect,state change,binding/association,dissociation,missing interaction target_generations: 3 to the web application of the project followed by KGML file hsa05200.xml.

### output examples

#### gravis
![alt text](https://github.com/shockwave742/KEGG-pathway-gravis-visualizer/blob/main/examples/gravis_example.gif)

#### NetworkX
![alt text](https://github.com/shockwave742/KEGG-pathway-gravis-visualizer/blob/main/examples/NetworkX_example.png)


### contact
for any question or suggestions contact me via:
- email: amirreza.alise@gmail.com

