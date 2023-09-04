from django.db import models

class ProteinSequence(models.Model):
    sequence = models.TextField()
    # Add any other fields you need for analysis results

    def __str__(self):
        return f"Protein Sequence {self.id}"

    class Meta:
        app_label = 'protein_app'