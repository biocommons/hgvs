import unittest

import hgvs.edit
from hgvs.exceptions import HGVSError

class Test_Edit(unittest.TestCase):
    def test_NARefAlt_exceptions(self):
        with self.assertRaises(HGVSError):
            edit = str(hgvs.edit.NARefAlt(None,None))

    def test_NARefAlt(self):
        self.assertEqual( str(hgvs.edit.NARefAlt('A','A'  ))  , '='             )
        self.assertEqual( str(hgvs.edit.NARefAlt('A','T'  ))  , 'A>T'           )
        self.assertEqual( str(hgvs.edit.NARefAlt('AA',None))  , 'delAA'         )
        self.assertEqual( str(hgvs.edit.NARefAlt(None,'TT'))  , 'insTT'         )
        self.assertEqual( str(hgvs.edit.NARefAlt('AA','T' ))  , 'delAAinsT'     )
        self.assertEqual( str(hgvs.edit.NARefAlt('A','TT' ))  , 'delAinsTT'     )


    def test_AARefAlt(self):
        self.assertEqual( str(hgvs.edit.AARefAlt('A','A'  ))  , '='             )
        self.assertEqual( str(hgvs.edit.AARefAlt('A','T'  ))  , 'Thr'           )
        self.assertEqual( str(hgvs.edit.AARefAlt('AA',None))  , 'del'           )
        self.assertEqual( str(hgvs.edit.AARefAlt(None,'TT'))  , 'insThrThr'     )
        self.assertEqual( str(hgvs.edit.AARefAlt('','T' ))  ,   'delinsThr'     )
        self.assertEqual( str(hgvs.edit.AARefAlt('AA','T' ))  , 'delinsThr'     )
        self.assertEqual( str(hgvs.edit.AARefAlt('A','TT' ))  , 'delinsThrThr'  )

        self.assertEqual( str(hgvs.edit.AARefAlt('A','T','fs*6'))  , 'Thrfs*6'  )

    def test_AASub(self):
        self.assertEqual( str(hgvs.edit.AASub('A','T'  ))  , 'Thr'           )
        self.assertEqual( str(hgvs.edit.AASub('A','T','fs*6'))  , 'Thrfs*6'  )


    def test_Dup(self):
        self.assertEqual( str(hgvs.edit.Dup())				, 'dup' 		)
        self.assertEqual( str(hgvs.edit.Dup('T'))			, 'dupT' 		)
        
    def test_Repeat(self):
        self.assertEqual( str(hgvs.edit.Repeat('CAG',12,34)), 'CAG(12_34)' 	)

    def test_Repeat_exceptions(self):
        with self.assertRaises(HGVSError):
            edit = str(hgvs.edit.Repeat('CAG',34,12))

if __name__ == '__main__':
    unittest.main()

## <LICENSE>
## Copyright 2014 HGVS Contributors (https://bitbucket.org/invitae/hgvs)
## 
## Licensed under the Apache License, Version 2.0 (the "License");
## you may not use this file except in compliance with the License.
## You may obtain a copy of the License at
## 
##     http://www.apache.org/licenses/LICENSE-2.0
## 
## Unless required by applicable law or agreed to in writing, software
## distributed under the License is distributed on an "AS IS" BASIS,
## WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
## See the License for the specific language governing permissions and
## limitations under the License.
## </LICENSE>
