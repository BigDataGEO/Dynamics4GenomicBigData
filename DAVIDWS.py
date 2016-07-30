import sys

import logging
import traceback as tb
import suds.metrics as metrics
from suds import *
from suds.client import Client
from datetime import datetime


def getClient():
  url = 'https://david.ncifcrf.gov/webservice/services/DAVIDWebService?wsdl'
  client = Client(url)
  client.wsdl.services[0].setlocation('https://david.ncifcrf.gov/webservice/services/DAVIDWebService.DAVIDWebServiceHttpSoap11Endpoint/')
  client.service.authenticate('j142857z@uth.edu') 
  return client

def getSpecies(inputIds, idType, listName, listType):
  client = getClient()
  client.service.addList(inputIds, idType, listName, listType)
  return client.service.getSpecies()

def getChartReport(inputIds, idType, listName, listType, thd, count):
  client = getClient()
  client.service.addList(inputIds, idType, listName, listType)
  return client.service.getChartReport(thd, count)

def getTableReport(inputIds, idType, listName, listType):
  client = getClient()
  client.service.addList(inputIds, idType, listName, listType)
  return client.service.getTableReport()

def getConversionTypes():
  client = getClient()
  return client.service.getConversionTypes()




