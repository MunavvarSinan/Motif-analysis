from django.urls import path
from . import views

urlpatterns = [
path('analyze/', views.analyze_protein, name='analyze_protein'),
]