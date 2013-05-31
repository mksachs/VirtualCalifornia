#!/usr/bin/env python

import json
import datetime

class VCEncoder(json.JSONEncoder):
    def default(self, obj):
        """
        default method is used if there is an unexpected object type
        in our example obj argument will be Decimal('120.50') and datetime
        in this encoder we are converting all Decimal to float and datetime to str
        """
        if isinstance(obj, datetime.datetime):
            obj = str(obj)
        else:
            obj = super(VCEncoder, self).default(obj)
        return obj