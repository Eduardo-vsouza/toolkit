

class GTF(object):
    def __init__(self, ):


class GTFParser(object):
    def __init__(self, gtf):
        self.gtf = gtf

    def parse(self):
        with open(self.gtf, 'r') as handler:
            lines = handler.readlines()


