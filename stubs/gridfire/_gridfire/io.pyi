"""
GridFire io bindings
"""
from __future__ import annotations
__all__: list[str] = ['NetworkFileParser', 'ParsedNetworkData', 'SimpleReactionListFileParser']
class NetworkFileParser:
    pass
class ParsedNetworkData:
    pass
class SimpleReactionListFileParser(NetworkFileParser):
    def parse(self, filename: str) -> ParsedNetworkData:
        """
        Parse a simple reaction list file and return a ParsedNetworkData object.
        """
