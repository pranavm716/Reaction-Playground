<template>
  <PageHeader highlighted-word="Playground" other-words="Mode"/>
  <ChemDrawTool @smiles="storeSmiles"/>
  <v-btn :color="$route.meta.navBarColor" v-if="smiles" @click="displayMol">Run Reactions</v-btn>
  <v-img v-if="img_base64" :src="`data:image/png;base64, ${img_base64}`"></v-img>
</template>

<script>
import axios from 'axios'
import ChemDrawTool from '@/components/ChemDrawTool.vue'
import PageHeader from '@/components/PageHeader.vue'

export default {
  name: 'PlaygroundMode',
  components: { ChemDrawTool, PageHeader },
  data () {
    return {
      smiles: null,
      img_base64: null
    }
  },
  methods: {
    storeSmiles (smiles) {
      this.smiles = smiles
      this.img_base64 = null
    },
    async displayMol () {
      await axios.get('/mol/image', { params: { mol_smiles: this.smiles } })
        .then(res => {
          this.img_base64 = res.data
        })
        .catch(error => {
          console.log(error.response)
        })
    }
  }
}
</script>
