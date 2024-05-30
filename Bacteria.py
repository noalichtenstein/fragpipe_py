class Bacteria:
    def __init__(self, phylum, bacteria_class, order, family, genus, species):
        self.phylum = phylum
        self.bacteria_class = bacteria_class
        self.order = order
        self.family = family
        self.genus = genus
        self.species = species

    def get_phylum(self):
        return self.phylum

    def get_bacteria_class(self):
        return self.bacteria_class

    def get_order(self):
        return self.order

    def get_family(self):
        return self.family

    def get_genus(self):
        return self.genus

    def get_species(self):
        return self.species

    def set_phylum(self, value):
        self._phylum = value

    def set_bacteria_class(self, value):
        self._bacteria_class = value

    def set_order(self, value):
        self._order = value

    def set_family(self, value):
        self._family = value

    def set_genus(self, value):
        self._genus = value

    def set_species(self, value):
        self._species = value