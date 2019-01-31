#!/usr/bin/env python
from .MetaDisassembler import MetaDisassembler

def main():
    md = MetaDisassembler()
    md._get_args()
    md.input_query(md.query)
    md.disassemble()

    return True

if __name__ == "__main__":
    main()
