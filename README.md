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

### contact
for any question or suggestion you can contact me via:
- email: amirreza.alise@gmail.com

