Schema of the Response of Json Response from SierraPy
=====================================================

This folder contains the json response from the sierrapy. This will help to parse the data well. Each json data that contains child structures is subdivided into further json schema so that we can get the outputs in a better way. the json schemas will be put in this file for better understanding

Overall JSON schemas
--------------------

```json
{
  "$schema": "http://json-schema.org/draft-04/schema#",
  "type": "array",
  "items": {
    "type": "object",
    "properties": {
      "inputSequence": {
        "type": "object",
        "properties": {
          "header": {
            "type": "string"
          },
          "SHA512": {
            "¬type": "string"
          }
        }
      },
      "strain": {
        "type": "object",
        "properties": {
          "name": {
            "type": "string"
          }
        }
      },
      "subtypeText": {
        "type": "string"
      },
      "validationResults": {
        "type": "array",
        "items": {
          "type": "object",
          "properties": {
            "level": {
              "type": "string"
            },
            "message": {
              "type": "string"
            }
          },
          "required": [
            "level",
            "message"
          ]
        }
      },
      "alignedGeneSequences": {
        "type": "array",
        "items": {
          "type": "object",
          "properties": {
            "firstAA": {
              "type": "number"
            },
            "lastAA": {
              "type": "number"
            },
            "gene": {
              "type": "object",
              "properties": {
                "name": {
                  "type": "string"
                },
                "length": {
                  "type": "number"
                }
              }
            },
            "mutations": {
              "type": "array",
              "items": {
                "type": "object",
                "properties": {
                  "consensus": {
                    "type": "string"
                  },
                  "position": {
                    "type": "number"
                  },
                  "AAs": {
                    "type": "string"
                  },
                  "isInsertion": {
                    "type": "boolean"
                  },
                  "isDeletion": {
                    "type": "boolean"
                  },
                  "isApobecMutation": {
                    "type": "boolean"
                  },
                  "isApobecDRM": {
                    "type": "boolean"
                  },
                  "isUnusual": {
                    "type": "boolean"
                  },
                  "isSDRM": {
                    "type": "boolean"
                  },
                  "hasStop": {
                    "type": "boolean"
                  },
                  "primaryType": {
                    "type": "string"
                  },
                  "text": {
                    "type": "string"
                  }
                },
                "required": [
                  "consensus",
                  "position",
                  "AAs",
                  "isInsertion",
                  "isDeletion",
                  "isApobecMutation",
                  "isApobecDRM",
                  "isUnusual",
                  "isSDRM",
                  "hasStop",
                  "primaryType",
                  "text"
                ]
              }
            },
            "SDRMs": {
              "type": "array",
              "items": {
                "type": "object",
                "properties": {
                  "text": {
                    "type": "string"
                  }
                },
                "required": [
                  "text"
                ]
              }
            },
            "alignedNAs": {
              "type": "string"
            },
            "alignedAAs": {
              "type": "string"
            },
            "prettyPairwise": {
              "type": "object",
              "properties": {
                "positionLine": {
                  "type": "array",
                  "items": {
                    "type": "string"
                  }
                },
                "refAALine": {
                  "type": "array",
                  "items": {
                    "type": "string"
                  }
                },
                "alignedNAsLine": {
                  "type": "array",
                  "items": {
                    "type": "string"
                  }
                },
                "mutationLine": {
                  "type": "array",
                  "items": {
                    "type": "string"
                  }
                }
              }
            }
          },
          "required": [
            "firstAA",
            "lastAA",
            "gene",
            "mutations",
            "SDRMs",
            "alignedNAs",
            "alignedAAs",
            "prettyPairwise"
          ]
        }
      },
      "drugResistance": {
        "type": "array",
        "items": {
          "type": "object",
          "properties": {
            "version": {
              "type": "object",
              "properties": {
                "text": {
                  "type": "string"
                },
                "publishDate": {
                  "type": "string"
                }
              }
            },
            "gene": {
              "type": "object",
              "properties": {
                "name": {
                  "type": "string"
                }
              }
            },
            "drugScores": {
              "type": "array",
              "items": {
                "type": "object",
                "properties": {
                  "drugClass": {
                    "type": "object",
                    "properties": {
                      "name": {
                        "type": "string"
                      }
                    }
                  },
                  "drug": {
                    "type": "object",
                    "properties": {
                      "name": {
                        "type": "string"
                      },
                      "displayAbbr": {
                        "type": "string"
                      }
                    }
                  },
                  "score": {
                    "type": "number"
                  },
                  "level": {
                    "type": "number"
                  },
                  "partialScores": {
                    "type": "array",
                    "items": {
                      "type": "object",
                      "properties": {
                        "mutations": {
                          "type": "array",
                          "items": {
                            "type": "object",
                            "properties": {
                              "text": {
                                "type": "string"
                              },
                              "primaryType": {
                                "type": "string"
                              },
                              "comments": {
                                "type": "array",
                                "items": {
                                  "type": "object",
                                  "properties": {
                                    "triggeredAAs": {
                                      "type": "string"
                                    },
                                    "type": {
                                      "type": "string"
                                    },
                                    "text": {
                                      "type": "string"
                                    }
                                  },
                                  "required": [
                                    "triggeredAAs",
                                    "type",
                                    "text"
                                  ]
                                }
                              }
                            },
                            "required": [
                              "text",
                              "primaryType",
                              "comments"
                            ]
                          }
                        },
                        "score": {
                          "type": "number"
                        }
                      },
                      "required": [
                        "mutations",
                        "score"
                      ]
                    }
                  },
                  "text": {
                    "type": "string"
                  }
                },
                "required": [
                  "drugClass",
                  "drug",
                  "score",
                  "level",
                  "partialScores",
                  "text"
                ]
              }
            }
          },
          "required": [
            "version",
            "gene",
            "drugScores"
          ]
        }
      }
    },
    "required": [
      "inputSequence",
      "strain",
      "subtypeText",
      "validationResults",
      "alignedGeneSequences",
      "drugResistance"
    ]
  }
}
```
