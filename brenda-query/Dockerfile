FROM python:3
ADD requirements.txt / 
RUN pip install -r requirements.txt
ADD main.py /
ADD rxn_data.csv /
CMD [ "python", "./main.py" ]
