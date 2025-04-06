// SPDX-License-Identifier: MIT
pragma solidity ^0.8.19;

contract MoleculeRegistry {
    struct Molecule {
        uint256 index;
        string smiles;
    }

    mapping(uint256 => Molecule) public molecules;
    uint256 public moleculeCount;

    event MoleculeStored(uint256 indexed index, string smiles);

   
    function storeMolecule(uint256 index, string memory smiles) public {
        molecules[index] = Molecule(index, smiles);
        moleculeCount++;
        emit MoleculeStored(index, smiles);
    }

    
    function getMolecule(uint256 index) public view returns (string memory) {
        return molecules[index].smiles;
    }
}
